use std::{fs::read_to_string, path::Path};

use castep_cell_io::{CellDocument, CellParser, LatticeParam};

use chemrust_core::data::{
    atom::{Atoms, CoreAtomData},
    geom::coordinates::CoordData,
    lattice::{
        cell_param::{LatticeVectors, UnitCellParameters},
        CrystalModel, LatticeCell,
    },
};
use nalgebra::{Matrix3, Point3, Vector3};

use crate::{search_sites, SearchConfig, SearchReports, SiteIndex};

fn load_model(model_rel_path: &str) -> Option<LatticeCell> {
    let root_dir = env!("CARGO_MANIFEST_DIR");
    let cell_path = Path::new(root_dir).join(model_rel_path);
    let content = read_to_string(cell_path).unwrap();
    let cell_model: CellDocument = CellParser::from(content.as_str()).parse().unwrap();
    let lattice_param_block = cell_model.lattice();
    let mut indices = Vec::new();
    let mut symbols = Vec::new();
    let mut coordinates = Vec::new();
    let mut labels = Vec::new();
    if let LatticeParam::LatticeCart(lat_cart) = lattice_param_block.parameter() {
        let data = [
            lat_cart.a().to_vec(),
            lat_cart.b().to_vec(),
            lat_cart.c().to_vec(),
        ];
        let lattice_vec = LatticeVectors::new(Matrix3::from_vec(data.concat()));
        cell_model
            .ionic_positions()
            .positions()
            .iter()
            .enumerate()
            .for_each(|(id, pos)| {
                indices.push(id);
                symbols.push(pos.symbol());

                coordinates.push(CoordData::Fractional(Point3::from(pos.coordinate())));
                labels.push(None);
            });
        let atoms = Atoms::new(indices, symbols, coordinates, labels);
        let model = LatticeCell::new(lattice_vec, atoms);
        Some(model)
    } else {
        None
    }
}

#[test]
fn new_algo() {
    let model = load_model("../scanner_test_models/H2TP001.cell").unwrap();
    let cartesian_coords: Vec<Point3<f64>> = model
        .get_atom_data()
        .coords()
        .iter()
        .map(|cd| match cd {
            CoordData::Fractional(frac) => {
                model.get_cell_parameters().cell_tensor()
                    * frac.map(|v| {
                        if !(0.0..=1.0).contains(&v) {
                            v - v.floor()
                        } else {
                            v
                        }
                    })
            }
            CoordData::Cartesian(cart) => *cart,
        })
        .collect();
    let dist = 1.95164_f64;
    let site_index = SiteIndex::new(cartesian_coords);
    todo!()
}

#[test]
fn test_search() {
    let model = load_model("../scanner_test_models/H2TP001.cell").unwrap();
    let lattice_vec = model.get_cell_parameters();
    let points: Vec<Point3<f64>> = model
        .get_atom_data()
        .coords()
        .iter()
        .map(|cd| match cd {
            CoordData::Fractional(frac) => {
                lattice_vec.cell_tensor()
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
    let results: SearchReports = search_sites(&site_index, &search_config);
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
