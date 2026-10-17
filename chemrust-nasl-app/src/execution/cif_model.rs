//! Bridge between a `crystal-cif-io` 0.2 `CoreCif` data block and the
//! `chemrust-core` model traits the search job consumes.

use std::str::FromStr;

use castep_periodic_table::element::ElementSymbol;
use chemrust_core::data::{
    atom::CoreAtomData,
    geom::coordinates::CoordData,
    lattice::{CellConstants, CrystalModel, LatticeVectors, UnitCellParameters},
};
use crystal_cif_io::data_dict::core_cif::cell::CellAngle;
use crystal_cif_io::{
    AtomSite, Cell, CharString, CoreCif, DataBlockHeading, Number, Numeric, UnquotedString,
};
use nalgebra::Point3;
use thiserror::Error;

use chemrust_nasl::{CoordSite, Visualize};

/// A failure while turning a CIF data block into a search model.
#[derive(Debug, Error)]
pub enum CifModelError {
    #[error("the data block misses the {item} item")]
    MissingCellItem { item: &'static str },
    #[error("atom {index} carries no coordinates")]
    NoCoordinates { index: usize },
    #[error("atom {index} carries no element symbol")]
    MissingSymbol { index: usize },
    #[error("atom {index} carries the unknown element symbol {symbol}")]
    UnknownSymbol { index: usize, symbol: String },
}

/// Atom data read from the `atom_site` loop, in the `chemrust-core`
/// representation the search job consumes.
#[derive(Clone, Debug)]
struct CifAtomData {
    indices: Vec<usize>,
    symbols: Vec<ElementSymbol>,
    coords: Vec<CoordData>,
    labels: Vec<Option<String>>,
}

impl CoreAtomData for CifAtomData {
    fn indices_repr(&self) -> Vec<usize> {
        self.indices.clone()
    }

    fn symbols_repr(&self) -> Vec<ElementSymbol> {
        self.symbols.clone()
    }

    fn coords_repr(&self) -> Vec<CoordData> {
        self.coords.clone()
    }

    fn labels_repr(&self) -> Vec<Option<String>> {
        self.labels.clone()
    }
}

/// A CIF data block as a `chemrust-core` crystal model.
#[derive(Clone, Debug)]
pub struct CifModel {
    lattice: LatticeVectors,
    atoms: CifAtomData,
}

impl CifModel {
    /// Read the lattice and atom data out of one CIF data block.
    ///
    /// The cell items carry Angstrom lengths and degree angles. A
    /// missing angle defaults to 90. Each row needs fractional or
    /// Cartesian coordinates, and a known element symbol.
    pub fn from_block(block: &CoreCif) -> Result<Self, CifModelError> {
        Ok(Self {
            lattice: block_lattice(&block.cell)?,
            atoms: block_atoms(&block.atom_site)?,
        })
    }

    /// The lattice of the data block.
    pub fn lattice(&self) -> &LatticeVectors {
        &self.lattice
    }
}

impl CrystalModel for CifModel {
    fn get_cell_parameters(&self) -> &impl UnitCellParameters {
        &self.lattice
    }

    fn get_atom_data(&self) -> &impl CoreAtomData {
        &self.atoms
    }

    fn get_cell_parameters_mut(&mut self) -> &mut impl UnitCellParameters {
        &mut self.lattice
    }

    fn get_atom_data_mut(&mut self) -> &mut impl CoreAtomData {
        &mut self.atoms
    }
}

/// Build the export CIF data block from a model and one new site.
///
/// The `cell` items come from the model lattice. The `atom_site`
/// rows come from the model atoms plus the new adsorbate site. The
/// rows carry labels, symbols, and coordinates only.
pub fn export_block(
    base_model: &impl CrystalModel,
    coord_site: &(impl CoordSite + Visualize),
    element_symbol: &ElementSymbol,
    block_name: &str,
) -> CoreCif {
    let cell = base_model.get_cell_parameters();
    let atom_data = base_model.get_atom_data();
    let tensor = cell.lattice_bases();
    let constants: CellConstants = CellConstants::from(&tensor);

    let cell_items = Cell {
        length_a: Some(Numeric::from(constants.length_a())),
        length_b: Some(Numeric::from(constants.length_b())),
        length_c: Some(Numeric::from(constants.length_c())),
        angle_alpha: Some(constants.angle_alpha().to_degrees().into()),
        angle_beta: Some(constants.angle_beta().to_degrees().into()),
        angle_gamma: Some(constants.angle_gamma().to_degrees().into()),
        ..Cell::default()
    };

    let new_index = atom_data.indices_repr().len();
    let mut rows: Vec<AtomSite> = atom_data
        .symbols_repr()
        .iter()
        .zip(atom_data.coords_repr().iter())
        .zip(atom_data.labels_repr().iter())
        .map(|((symbol, coord), label)| {
            let base = AtomSite {
                type_symbol: Some(char_string(symbol.to_string())),
                label: label.as_deref().map(|text| char_string(text.to_string())),
                ..AtomSite::default()
            };
            match coord {
                CoordData::Fractional(frac) => AtomSite {
                    fract_x: Some(Numeric::from(frac.x)),
                    fract_y: Some(Numeric::from(frac.y)),
                    fract_z: Some(Numeric::from(frac.z)),
                    ..base
                },
                CoordData::Cartesian(cart) => AtomSite {
                    cartn_x: Some(Numeric::from(cart.x)),
                    cartn_y: Some(Numeric::from(cart.y)),
                    cartn_z: Some(Numeric::from(cart.z)),
                    ..base
                },
            }
        })
        .collect();

    let frac = coord_site.fractional_coord(tensor);
    let new_row = AtomSite {
        type_symbol: Some(char_string(element_symbol.to_string())),
        label: Some(char_string(format!("{}{}", element_symbol, new_index))),
        fract_x: Some(Numeric::from(frac.x)),
        fract_y: Some(Numeric::from(frac.y)),
        fract_z: Some(Numeric::from(frac.z)),
        ..AtomSite::default()
    };
    rows.push(new_row);

    CoreCif {
        heading: DataBlockHeading::new(block_name.to_string()),
        cell: cell_items,
        atom_site: rows,
        ..CoreCif::default()
    }
}

fn block_lattice(cell: &Cell) -> Result<LatticeVectors, CifModelError> {
    let a = cell_length(cell, "cell_length_a")?;
    let b = cell_length(cell, "cell_length_b")?;
    let c = cell_length(cell, "cell_length_c")?;
    let constants = CellConstants::new(
        a,
        b,
        c,
        cell_angle(cell.angle_alpha.as_ref()),
        cell_angle(cell.angle_beta.as_ref()),
        cell_angle(cell.angle_gamma.as_ref()),
    );
    Ok(LatticeVectors::new(constants.lattice_bases()))
}

fn cell_length(cell: &Cell, item: &'static str) -> Result<f64, CifModelError> {
    let value = match item {
        "cell_length_a" => &cell.length_a,
        "cell_length_b" => &cell.length_b,
        _ => &cell.length_c,
    };
    value
        .as_ref()
        .map(numeric_f64)
        .ok_or(CifModelError::MissingCellItem { item })
}

fn cell_angle(angle: Option<&CellAngle>) -> f64 {
    angle.map(|v| v.angle().0).unwrap_or(90.0).to_radians()
}

fn block_atoms(rows: &[AtomSite]) -> Result<CifAtomData, CifModelError> {
    let mut symbols: Vec<ElementSymbol> = Vec::with_capacity(rows.len());
    let mut coords: Vec<CoordData> = Vec::with_capacity(rows.len());
    let mut labels: Vec<Option<String>> = Vec::with_capacity(rows.len());
    for (index, row) in rows.iter().enumerate() {
        let text = row
            .type_symbol
            .as_ref()
            .map(|symbol| symbol.content())
            .ok_or(CifModelError::MissingSymbol { index })?;
        let symbol = ElementSymbol::from_str(text).map_err(|_| CifModelError::UnknownSymbol {
            index,
            symbol: text.to_string(),
        })?;
        symbols.push(symbol);
        coords.push(row_coords(row, index)?);
        labels.push(row.label.as_ref().map(|text| text.content().to_string()));
    }
    Ok(CifAtomData {
        indices: (0..rows.len()).collect(),
        symbols,
        coords,
        labels,
    })
}

fn row_coords(row: &AtomSite, index: usize) -> Result<CoordData, CifModelError> {
    if let (Some(x), Some(y), Some(z)) = (&row.fract_x, &row.fract_y, &row.fract_z) {
        return Ok(CoordData::Fractional(Point3::new(
            numeric_f64(x),
            numeric_f64(y),
            numeric_f64(z),
        )));
    }
    if let (Some(x), Some(y), Some(z)) = (&row.cartn_x, &row.cartn_y, &row.cartn_z) {
        return Ok(CoordData::Cartesian(Point3::new(
            numeric_f64(x),
            numeric_f64(y),
            numeric_f64(z),
        )));
    }
    Err(CifModelError::NoCoordinates { index })
}

fn numeric_f64(value: &Numeric) -> f64 {
    match value.number() {
        Number::Float(v) => v.0,
        Number::Integer(v) => v.0 as f64,
    }
}

fn char_string(text: String) -> CharString {
    CharString::Unquoted(UnquotedString::new(text))
}
