use std::{fs::read_to_string, path::Path};

use castep_cell_io::{CellDocument, CellParseError, CellParser};

use chemrust_core::data::atom::CoreAtomData;
use chemrust_core::data::geom::coordinates::CoordData;
use chemrust_core::data::lattice::CrystalModel;
use chemrust_core::data::lattice::UnitCellParameters;
use kd_tree::KdIndexTree;
use nalgebra::{Point3, Vector3};

use crate::{search_sites, SearchConfig, SearchReports, SiteIndex};

fn load_model(model_rel_path: &str) -> Result<CellDocument, CellParseError> {
    let root_dir = env!("CARGO_MANIFEST_DIR");
    let cell_path = Path::new(root_dir).join(model_rel_path);
    let content = read_to_string(cell_path).unwrap();
    CellParser::from(&content).parse()
}

#[test]
fn test_search() {
    let model = load_model("../scanner_test_models/g-C3N4.cell").unwrap();
    let lattice_vec = model.get_cell_parameters();
    let points: Vec<Point3<f64>> = model
        .get_atom_data()
        .coords_repr()
        .iter()
        .map(|cd| match cd {
            CoordData::Fractional(frac) => {
                lattice_vec.lattice_bases()
                    * frac.map(|v| {
                        if !(0.0..=1.0).contains(&v) {
                            v - v.floor()
                        } else {
                            v
                        }
                    })
            }
            CoordData::Cartesian(p) => *p,
        })
        .collect();
    let dist: f64 = 1.95164;
    let site_index = SiteIndex::new(&points);
    let search_points: Vec<(usize, Point3<f64>)> =
        points.iter().enumerate().map(|(i, p)| (i, *p)).collect();
    let search_config = SearchConfig::new(&search_points, dist);
    let _results: SearchReports = search_sites(&site_index, &search_config);
}
