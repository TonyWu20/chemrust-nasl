//! Converters between the `chemrust-geometry` `Structure` model and the
//! two file formats the app reads and writes: CASTEP `.cell` 0.7 and
//! `crystal-cif-io` 0.2 data blocks.

use std::str::FromStr;

use castep_cell_fmt::{format::to_string_many_spaced, ToCellFile};
use castep_cell_io::{
    cell::{
        lattice_param::LatticeCart,
        positions::{PositionFracEntry, PositionsFrac},
        species::Species,
    },
    CellDocument, Lattice, Positions,
};
use castep_periodic_table::element::ElementSymbol;
use chemrust_geometry::{CellConstants, FracCoord, LatticeVectors, Structure};
use crystal_cif_io::data_dict::core_cif::cell::CellAngle;
use crystal_cif_io::{
    AtomSite, Cell, CharString, CoreCif, DataBlockHeading, Number, Numeric, UnquotedString,
};
use nalgebra::Matrix3;
use thiserror::Error;

use chemrust_nasl::{CoordSite, Visualize};

/// A failure while turning a file model into a `Structure`.
#[derive(Debug, Error)]
pub enum ModelError {
    #[error("the data block misses the {item} item")]
    MissingCellItem { item: &'static str },
    #[error("atom {index} carries no coordinates")]
    NoCoordinates { index: usize },
    #[error("atom {index} carries no element symbol")]
    MissingSymbol { index: usize },
    #[error("atom {index} carries the unknown element symbol {symbol}")]
    UnknownSymbol { index: usize, symbol: String },
}

fn lattice_from_cell_items(cell: &Cell) -> Result<LatticeVectors, ModelError> {
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
    Ok(LatticeVectors::from_constants(constants))
}

fn cell_length(cell: &Cell, item: &'static str) -> Result<f64, ModelError> {
    let value = match item {
        "cell_length_a" => &cell.length_a,
        "cell_length_b" => &cell.length_b,
        _ => &cell.length_c,
    };
    value
        .as_ref()
        .map(numeric_f64)
        .ok_or(ModelError::MissingCellItem { item })
}

fn cell_angle(angle: Option<&CellAngle>) -> f64 {
    angle.map(|v| v.angle().0).unwrap_or(90.0).to_radians()
}

/// Read a `Structure` out of a CASTEP `.cell` document.
///
/// `LATTICE_CART` rows are the a, b, c vectors; `LATTICE_ABC` angles are
/// in degrees. Position entries carry a species identifier, read as a
/// castep `ElementSymbol`.
pub fn cell_to_structure(doc: &CellDocument) -> Result<Structure, ModelError> {
    let lattice = match &doc.lattice {
        Lattice::Cart(cart) => LatticeVectors::new(Matrix3::from_columns(&[
            cart.a.into(),
            cart.b.into(),
            cart.c.into(),
        ])),
        Lattice::Abc(abc) => LatticeVectors::from_constants(CellConstants::new(
            abc.abc[0],
            abc.abc[1],
            abc.abc[2],
            abc.angles[0].to_radians(),
            abc.angles[1].to_radians(),
            abc.angles[2].to_radians(),
        )),
    };
    let tensor = lattice.tensor();
    let mut species: Vec<ElementSymbol> = Vec::new();
    let mut frac_coords: Vec<FracCoord> = Vec::new();
    match &doc.positions {
        Positions::Frac(frac) => {
            for (index, entry) in frac.positions.iter().enumerate() {
                species.push(species_symbol(&entry.species, index)?);
                frac_coords.push(FracCoord::from_array(entry.coord));
            }
        }
        Positions::Abs(abs) => {
            let inverse = tensor.try_inverse().ok_or(ModelError::MissingCellItem {
                item: "LATTICE_CART",
            })?;
            for (index, entry) in abs.positions.iter().enumerate() {
                species.push(species_symbol(&entry.species, index)?);
                let cart = nalgebra::Point3::new(entry.coord[0], entry.coord[1], entry.coord[2]);
                let frac = inverse * cart;
                frac_coords.push(FracCoord::new(frac.x, frac.y, frac.z));
            }
        }
    }
    let count = species.len();
    Ok(Structure::new(
        species,
        frac_coords,
        Some(lattice),
        [true; 3],
        vec![0; count],
        vec![None; count],
        None,
    ))
}

fn species_symbol(sp: &Species, index: usize) -> Result<ElementSymbol, ModelError> {
    if let Some(symbol) = sp.as_symbol() {
        return ElementSymbol::from_str(symbol).map_err(|_| ModelError::UnknownSymbol {
            index,
            symbol: symbol.clone(),
        });
    }
    let z = sp
        .as_atomic_number()
        .ok_or(ModelError::MissingSymbol { index })?;
    ElementSymbol::try_from(*z).map_err(|_| ModelError::UnknownSymbol {
        index,
        symbol: z.to_string(),
    })
}

/// Read a `Structure` out of one `crystal-cif-io` 0.2 data block.
///
/// Cell items carry Angstrom lengths and degree angles. A missing angle
/// defaults to 90. Each row needs a known element symbol.
pub fn cif_to_structure(block: &CoreCif) -> Result<Structure, ModelError> {
    let lattice = lattice_from_cell_items(&block.cell)?;
    let tensor = lattice.tensor();
    let inverse = tensor.try_inverse().ok_or(ModelError::MissingCellItem {
        item: "LATTICE_CART",
    })?;
    let mut species: Vec<ElementSymbol> = Vec::with_capacity(block.atom_site.len());
    let mut frac_coords: Vec<FracCoord> = Vec::with_capacity(block.atom_site.len());
    let mut labels: Vec<Option<String>> = Vec::with_capacity(block.atom_site.len());
    for (index, row) in block.atom_site.iter().enumerate() {
        let text = row
            .type_symbol
            .as_ref()
            .map(|symbol| symbol.content())
            .ok_or(ModelError::MissingSymbol { index })?;
        species.push(
            ElementSymbol::from_str(text).map_err(|_| ModelError::UnknownSymbol {
                index,
                symbol: text.to_string(),
            })?,
        );
        frac_coords.push(row_coords(row, index, &inverse)?);
        labels.push(row.label.as_ref().map(|text| text.content().to_string()));
    }
    let count = species.len();
    Ok(Structure::new(
        species,
        frac_coords,
        Some(lattice),
        [true; 3],
        vec![0; count],
        labels,
        None,
    ))
}

fn row_coords(
    row: &AtomSite,
    index: usize,
    inverse: &Matrix3<f64>,
) -> Result<FracCoord, ModelError> {
    if let (Some(x), Some(y), Some(z)) = (&row.fract_x, &row.fract_y, &row.fract_z) {
        return Ok(FracCoord::new(
            numeric_f64(x),
            numeric_f64(y),
            numeric_f64(z),
        ));
    }
    if let (Some(x), Some(y), Some(z)) = (&row.cartn_x, &row.cartn_y, &row.cartn_z) {
        let cart = nalgebra::Point3::new(numeric_f64(x), numeric_f64(y), numeric_f64(z));
        let frac = inverse * cart;
        return Ok(FracCoord::new(frac.x, frac.y, frac.z));
    }
    Err(ModelError::NoCoordinates { index })
}

/// Build the export CIF text from a model and one new adsorbate site.
///
/// The `cell` items come from the model lattice. The `atom_site` rows
/// come from the model atoms plus the new site. Rows carry labels,
/// symbols, and coordinates only.
pub fn cif_export_text(
    base_model: &Structure,
    coord_site: &(impl CoordSite + Visualize),
    element_symbol: &ElementSymbol,
    block_name: &str,
) -> String {
    let tensor = model_tensor(base_model);
    let block = cif_export_block(base_model, coord_site, element_symbol, block_name, tensor);
    crystal_cif_io::write_core_cif(std::slice::from_ref(&block))
}

/// Build the export `.cell` text from a model and one new site.
pub fn cell_export_text(
    base_model: &Structure,
    coord_site: &(impl CoordSite + Visualize),
    element_symbol: &ElementSymbol,
) -> String {
    let tensor = model_tensor(base_model);
    let doc = cell_export_document(base_model, coord_site, element_symbol, tensor)
        .expect("the structure converts to a cell document");
    to_string_many_spaced(&doc.to_cell_file())
}

fn cif_export_block(
    base_model: &Structure,
    coord_site: &(impl CoordSite + Visualize),
    element_symbol: &ElementSymbol,
    block_name: &str,
    tensor: &Matrix3<f64>,
) -> CoreCif {
    let cell_vectors = model_cell(base_model);
    let constants = CellConstants::from(cell_vectors);
    let cell_items = Cell {
        length_a: Some(Numeric::from(constants.a)),
        length_b: Some(Numeric::from(constants.b)),
        length_c: Some(Numeric::from(constants.c)),
        angle_alpha: Some(constants.alpha.to_degrees().into()),
        angle_beta: Some(constants.beta.to_degrees().into()),
        angle_gamma: Some(constants.gamma.to_degrees().into()),
        ..Cell::default()
    };

    let new_index = base_model.species.len();
    let rows: Vec<AtomSite> = base_model
        .species
        .iter()
        .zip(base_model.frac_coords.iter())
        .zip(base_model.labels.iter())
        .map(|((symbol, coord), label)| {
            let base = AtomSite {
                type_symbol: Some(char_string(symbol.to_string())),
                label: label.as_deref().map(|text| char_string(text.to_string())),
                ..AtomSite::default()
            };
            AtomSite {
                fract_x: Some(Numeric::from(coord.0.x)),
                fract_y: Some(Numeric::from(coord.0.y)),
                fract_z: Some(Numeric::from(coord.0.z)),
                ..base
            }
        })
        .collect();

    let frac = coord_site.fractional_coord(*tensor);
    let new_row = AtomSite {
        type_symbol: Some(char_string(element_symbol.to_string())),
        label: Some(char_string(format!("{}{}", element_symbol, new_index))),
        fract_x: Some(Numeric::from(frac.x)),
        fract_y: Some(Numeric::from(frac.y)),
        fract_z: Some(Numeric::from(frac.z)),
        ..AtomSite::default()
    };
    let mut rows = rows;
    rows.push(new_row);

    CoreCif {
        heading: DataBlockHeading::new(block_name.to_string()),
        cell: cell_items,
        atom_site: rows,
        ..CoreCif::default()
    }
}

fn cell_export_document(
    base_model: &Structure,
    coord_site: &(impl CoordSite + Visualize),
    element_symbol: &ElementSymbol,
    tensor: &Matrix3<f64>,
) -> Result<CellDocument, castep_cell_fmt::Error> {
    let lattice = Lattice::Cart(LatticeCart {
        unit: None,
        a: tensor.column(0).into(),
        b: tensor.column(1).into(),
        c: tensor.column(2).into(),
    });
    let entries: Vec<PositionFracEntry> = base_model
        .species
        .iter()
        .zip(base_model.frac_coords.iter())
        .map(|(symbol, coord)| {
            PositionFracEntry::builder()
                .species(Species::Symbol(symbol.to_string()))
                .coord(coord.into_array())
                .build()
        })
        .collect();
    let frac = coord_site.fractional_coord(*tensor);
    let mut entries = entries;
    entries.push(
        PositionFracEntry::builder()
            .species(Species::Symbol(element_symbol.to_string()))
            .coord([frac.x, frac.y, frac.z])
            .build(),
    );
    CellDocument::builder()
        .lattice(lattice)
        .positions(Positions::Frac(PositionsFrac { positions: entries }))
        .build()
}

/// The lattice of a periodic model.
fn model_cell(model: &Structure) -> LatticeVectors {
    model.cell.expect("the model carries a lattice")
}

/// The lattice tensor of a periodic model.
fn model_tensor(model: &Structure) -> &Matrix3<f64> {
    model
        .cell
        .as_ref()
        .expect("the model carries a lattice")
        .tensor()
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
