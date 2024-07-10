use std::{fs::read_to_string, path::Path};

use castep_cell_io::{CellDocument, CellParser};

use chemrust_core::data::{
    atom::CoreAtomData,
    geom::coordinates::CoordData,
    lattice::{CrystalModel, UnitCellParameters},
};
use nalgebra::{Point3, Vector3};

use crate::{AdsSiteLocator, SearchConfig, SearchReports, SiteIndex};

fn load_model(model_rel_path: &str) -> Option<CellDocument> {
    let root_dir = env!("CARGO_MANIFEST_DIR");
    let cell_path = Path::new(root_dir).join(model_rel_path);
    let content = read_to_string(cell_path).unwrap();
    CellParser::from(&content).parse().ok()
}

#[test]
fn test_search() {
    let model = load_model("../scanner_test_models/H2TP001.cell").unwrap();
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
    let site_index = SiteIndex::new(points);
    let search_points: Vec<(usize, Point3<f64>)> = site_index
        .coords()
        .iter()
        .enumerate()
        .map(|(i, p)| (i, *p))
        .collect();
    let search_config = SearchConfig::new(&search_points, dist);
    let searcher = AdsSiteLocator::new(&site_index, &search_config);
    let results: SearchReports = searcher.search_sites();
    results
        .points()
        .unwrap()
        .iter()
        .find(|p| p.atom_ids() == [10, 12, 29, 31].to_vec())
        .map(|special| {
            special.atom_ids().iter().for_each(|id| {
                let p = site_index.coords()[*id];
                let od: Vector3<f64> = special.point() - p;
                println!("Distance with {id}: {}", od.norm());
            });
        })
        .expect("No center");
}
