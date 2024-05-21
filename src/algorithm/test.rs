use std::{fs::read_to_string, path::Path};

use castep_cell_io::{CellDocument, CellParser, IonicPosition, IonicPositionBlock, LatticeParam};

use castep_periodic_table::element::ElementSymbol;
use chemrust_core::data::{
    atom::{Atoms, CoreAtomData},
    geom::coordinates::CoordData,
    lattice::{cell_param::LatticeVectors, CrystalModel, LatticeCell},
};
use nalgebra::{Matrix3, Point3};

use crate::coordination_sites::Visualize;

use super::{search_sites, SearchConfig, SiteIndex};

#[test]
fn test_search() {
    let root_dir = env!("CARGO_MANIFEST_DIR");
    let cell_path = Path::new(root_dir).join("SAC_GDY_V.cell");
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
        let points: Vec<Point3<f64>> = model
            .get_atom_data()
            .coords()
            .iter()
            .map(|cd| match cd {
                CoordData::Fractional(frac) => {
                    let coord: Vec<f64> = frac
                        .iter()
                        .map(|&v| {
                            if v < 0.0 {
                                v + 1.0
                            } else if v > 1.0 {
                                v - 1.0
                            } else {
                                v
                            }
                        })
                        .collect();
                    lattice_vec.tensor() * Point3::from_slice(&coord)
                }
                CoordData::Cartesian(p) => *p,
            })
            .collect();
        let dist: f64 = 2.6;
        let site_index = SiteIndex::new(points);
        let search_points: Vec<(usize, Point3<f64>)> = site_index
            .coords()
            .iter()
            .enumerate()
            .map(|(i, p)| (i, *p))
            .collect();
        let search_config = SearchConfig::new(&search_points, dist);
        let results = search_sites(&site_index, &search_config);
        println!(
            "Spheres: {}, Circles: {}, Points: {}",
            results.spheres.len(),
            results.circles.len(),
            results.points.len()
        );
        results.circles.iter().for_each(|coord_circle| {
            let atom_ids_text = coord_circle
                .atom_ids
                .iter()
                .map(|id| id.to_string())
                .collect::<Vec<String>>()
                .join("_");
            let atom = coord_circle.draw_with_element(&ElementSymbol::Pt);
            let coordinate = lattice_vec.tensor().try_inverse().unwrap() * atom.coord().xyz();
            let new_point = IonicPosition::new(ElementSymbol::Pt, coordinate.into(), None);
            let mut new_positions = cell_model.ionic_positions().positions().to_vec();
            new_positions.push(new_point);
            let new_positions_block = IonicPositionBlock::new(
                cell_model.ionic_positions().unit(),
                new_positions,
                cell_model.ionic_positions().keyword(),
                cell_model.ionic_positions().spin_polarised(),
            );
            let new_model = CellDocument::new(lattice_param_block, new_positions_block);
            let filename = format!("demo/double/SAC_GDY_V_Pt_{dist}_{}.cell", atom_ids_text);
            new_model.write_out(filename).expect("Write out error")
        });
        results.points.iter().for_each(|coord_point| {
            let atom_ids_text = coord_point
                .atom_ids
                .iter()
                .map(|id| id.to_string())
                .collect::<Vec<String>>()
                .join("_");
            let coordinate = lattice_vec.tensor().try_inverse().unwrap() * coord_point.point;
            let new_point = IonicPosition::new(ElementSymbol::Pt, coordinate.into(), None);
            let mut new_positions = cell_model.ionic_positions().positions().to_vec();
            new_positions.push(new_point);
            let new_positions_block = IonicPositionBlock::new(
                cell_model.ionic_positions().unit(),
                new_positions,
                cell_model.ionic_positions().keyword(),
                cell_model.ionic_positions().spin_polarised(),
            );
            let new_model = CellDocument::new(lattice_param_block, new_positions_block);
            let filename = format!("demo/multi/SAC_GDY_V_Pt_{dist}_{}.cell", atom_ids_text);
            new_model.write_out(filename).expect("Write out error")
        })
    }
}
