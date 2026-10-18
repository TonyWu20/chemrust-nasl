#![allow(dead_code)]
//! Model loading shared by the benchmark examples.
//!
//! `.cell` files come from `castep-cell-fmt` 0.3 plus `castep-cell-io`
//! 0.7. `.cif` files come from `crystal-cif-io` 0.2. Both paths build
//! the `chemrust-geometry` `Structure` model.

use std::fs::read_to_string;
use std::path::Path;
use std::str::FromStr;

use castep_cell_fmt::parse as parse_cell;
use castep_cell_io::{
    cell::species::Species, CellDocument, Lattice, Positions,
};
use chemrust_geometry::{CellConstants, FracCoord, LatticeVectors, Structure};
use crystal_cif_io::{parse_core_cif, Number, Numeric};
use nalgebra::{Matrix3, Point3, Vector3};

/// Load a `.cell` or `.cif` model into one periodic `Structure`.
pub fn load_structure(path: &Path) -> Result<Structure, String> {
    let text = read_to_string(path).map_err(|e| e.to_string())?;
    let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("");
    match ext {
        "cell" => cell_to_structure(&parse_cell(&text).map_err(|e| e.to_string())?),
        "cif" => cif_to_structure(&text),
        _ => Err(format!("unsupported model file: {path:?}")),
    }
}

/// Cartesian points (Angstrom) of every atom, wrapped into the cell.
pub fn cart_points(s: &Structure) -> Vec<Point3<f64>> {
    let tensor = s
        .cell
        .as_ref()
        .expect("the model carries a lattice")
        .tensor();
    s.frac_coords
        .iter()
        .map(|c| {
            let mut f = *c;
            f.wrap();
            tensor * Point3::new(f.0.x, f.0.y, f.0.z)
        })
        .collect()
}

fn cell_to_structure(doc: &CellDocument) -> Result<Structure, String> {
    let lattice = match &doc.lattice {
        Lattice::Cart(cart) => LatticeVectors::new(Matrix3::from_columns(&[
            Vector3::from(cart.a),
            Vector3::from(cart.b),
            Vector3::from(cart.c),
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
    let inverse = tensor.try_inverse().ok_or("singular lattice")?;
    let mut species = Vec::new();
    let mut coords = Vec::new();
    match &doc.positions {
        Positions::Frac(frac) => {
            for entry in frac.positions.iter() {
                species.push(cell_species(&entry.species)?);
                coords.push(FracCoord::from_array(entry.coord));
            }
        }
        Positions::Abs(abs) => {
            for entry in abs.positions.iter() {
                species.push(cell_species(&entry.species)?);
                let cart = Point3::new(entry.coord[0], entry.coord[1], entry.coord[2]);
                let frac = inverse * cart;
                coords.push(FracCoord::new(frac.x, frac.y, frac.z));
            }
        }
    }
    let n = species.len();
    Ok(Structure::new(
        species,
        coords,
        Some(lattice),
        [true; 3],
        vec![0; n],
        vec![None; n],
        None,
    ))
}

fn cell_species(sp: &Species) -> Result<castep_periodic_table::element::ElementSymbol, String> {
    use castep_periodic_table::element::ElementSymbol;
    if let Some(symbol) = sp.as_symbol() {
        return ElementSymbol::from_str(symbol)
            .map_err(|_| format!("unknown element symbol {symbol}"));
    }
    let z = sp
        .as_atomic_number()
        .ok_or("position entry carries no species")?;
    ElementSymbol::try_from(*z).map_err(|_| format!("unknown atomic number {z}"))
}

fn cif_to_structure(text: &str) -> Result<Structure, String> {
    let doc = parse_core_cif(text).map_err(|e| e.to_string())?;
    let block = doc.blocks.first().cloned().ok_or("no data block")?;
    let cell = &block.cell;
    let a = cif_len(&cell.length_a, "cell_length_a")?;
    let b = cif_len(&cell.length_b, "cell_length_b")?;
    let c = cif_len(&cell.length_c, "cell_length_c")?;
    let angle = |v: &Option<crystal_cif_io::data_dict::core_cif::cell::CellAngle>, default: f64| {
        v.as_ref()
            .map(|v| v.angle().0)
            .unwrap_or(default)
            .to_radians()
    };
    let lattice = LatticeVectors::from_constants(CellConstants::new(
        a,
        b,
        c,
        angle(&cell.angle_alpha, 90.0),
        angle(&cell.angle_beta, 90.0),
        angle(&cell.angle_gamma, 90.0),
    ));
    let tensor = lattice.tensor();
    let inverse = tensor.try_inverse().ok_or("singular lattice")?;
    let mut species = Vec::with_capacity(block.atom_site.len());
    let mut coords = Vec::with_capacity(block.atom_site.len());
    for (i, row) in block.atom_site.iter().enumerate() {
        use castep_periodic_table::element::ElementSymbol;
        let symbol = row
            .type_symbol
            .as_ref()
            .map(|s| s.content())
            .ok_or_else(|| format!("atom {i} carries no element symbol"))?;
        species.push(
            ElementSymbol::from_str(symbol)
                .map_err(|_| format!("unknown element symbol {symbol}"))?,
        );
        if let (Some(x), Some(y), Some(z)) = (&row.fract_x, &row.fract_y, &row.fract_z) {
            coords.push(FracCoord::new(
                numeric_f64(x),
                numeric_f64(y),
                numeric_f64(z),
            ));
        } else if let (Some(x), Some(y), Some(z)) =
            (&row.cartn_x, &row.cartn_y, &row.cartn_z)
        {
            let cart = Point3::new(numeric_f64(x), numeric_f64(y), numeric_f64(z));
            let frac = inverse * cart;
            coords.push(FracCoord::new(frac.x, frac.y, frac.z));
        } else {
            return Err(format!("atom {i} carries no coordinates"));
        }
    }
    let n = species.len();
    Ok(Structure::new(
        species,
        coords,
        Some(lattice),
        [true; 3],
        vec![0; n],
        vec![None; n],
        None,
    ))
}

fn cif_len(v: &Option<Numeric>, item: &str) -> Result<f64, String> {
    v.as_ref()
        .map(numeric_f64)
        .ok_or_else(|| format!("missing {item}"))
}

fn numeric_f64(value: &Numeric) -> f64 {
    match value.number() {
        Number::Float(v) => v.0,
        Number::Integer(v) => v.0 as f64,
    }
}
