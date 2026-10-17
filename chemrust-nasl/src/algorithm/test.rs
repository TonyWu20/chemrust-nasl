use std::{fs::read_to_string, path::Path};

use castep_cell_fmt::parse;
use castep_cell_io::{CellDocument, Lattice, Positions};
use chemrust_geometry::{CellConstants, LatticeVectors};
use nalgebra::{Matrix3, Point3, Vector3};

use crate::{SearchConfig, SearchReports};

/// Wrap one fractional coordinate into `[0.0, 1.0)` the same way the
/// legacy test did.
fn wrap_frac(v: f64) -> f64 {
    if !(0.0..=1.0).contains(&v) {
        v - v.floor()
    } else {
        v
    }
}

fn load_model(model_rel_path: &str) -> Result<CellDocument, castep_cell_fmt::Error> {
    let root_dir = env!("CARGO_MANIFEST_DIR");
    let cell_path = Path::new(root_dir).join(model_rel_path);
    let content = read_to_string(cell_path).unwrap();
    parse(&content)
}

/// Lattice vectors in the `chemrust-geometry` convention: the columns
/// of the tensor are the a, b and c vectors in Cartesian coordinates.
fn lattice_vectors(doc: &CellDocument) -> LatticeVectors {
    match &doc.lattice {
        // LATTICE_CART rows are the a, b and c vectors.
        Lattice::Cart(cart) => LatticeVectors::new(Matrix3::from_columns(&[
            Vector3::from(cart.a),
            Vector3::from(cart.b),
            Vector3::from(cart.c),
        ])),
        // LATTICE_ABC angles are in degrees; `CellConstants` uses radians.
        Lattice::Abc(abc) => LatticeVectors::from_constants(CellConstants::new(
            abc.abc[0],
            abc.abc[1],
            abc.abc[2],
            abc.angles[0].to_radians(),
            abc.angles[1].to_radians(),
            abc.angles[2].to_radians(),
        )),
    }
}

/// Cartesian positions (Angstrom) of all atoms in the document.
fn cartesian_points(doc: &CellDocument) -> Vec<Point3<f64>> {
    let lattice = lattice_vectors(doc);
    let tensor = lattice.tensor();
    match &doc.positions {
        Positions::Frac(frac) => frac
            .positions
            .iter()
            .map(|entry| {
                let f = [
                    wrap_frac(entry.coord[0]),
                    wrap_frac(entry.coord[1]),
                    wrap_frac(entry.coord[2]),
                ];
                tensor * Point3::new(f[0], f[1], f[2])
            })
            .collect(),
        // Assumed Angstrom, the CASTEP default.
        Positions::Abs(abs) => abs
            .positions
            .iter()
            .map(|entry| Point3::new(entry.coord[0], entry.coord[1], entry.coord[2]))
            .collect(),
    }
}

#[test]
fn test_search() {
    let model = load_model("../scanner_test_models/SAC_GDY_V.cell").unwrap();
    let points = cartesian_points(&model);
    let dist: f64 = 1.95164;
    let search_points: Vec<(usize, Point3<f64>)> =
        points.iter().enumerate().map(|(i, p)| (i, *p)).collect();
    let search_config = SearchConfig::new(&search_points, &points, dist);
    let results: SearchReports = search_config.search_sites();
    // Site counts verified against the legacy `castep-cell-io` 0.2.18 +
    // `chemrust-core` 0.5.1 test: both stacks produce 73 / 93 / 26.
    assert_eq!(results.viable_single_points().unwrap().iter().len(), 73);
    assert_eq!(results.viable_double_points().unwrap().iter().len(), 93);
    assert_eq!(results.points().unwrap().iter().len(), 26);
}
