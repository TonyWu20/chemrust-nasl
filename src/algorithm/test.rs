use std::{collections::HashSet, fs::read_to_string, path::Path};

use castep_cell_parser::{CellDocument, CellParser, LatticeParam};

use castep_periodic_table::element::ElementSymbol;
use chemrust_core::data::{
    atom::CoreAtomData,
    geom::coordinates::CoordData,
    lattice::{
        cell_param::{LatticeVectors, UnitCellParameters},
        CrystalModel,
    },
};
use kiddo::SquaredEuclidean;
use nalgebra::{Matrix3, Point3};

use crate::{
    algorithm::{
        build_kd_tree_from_points,
        coordinate_sites::{CoordCircle, CoordPoint, CoordResult},
    },
    geometry::{Intersect, Sphere, SphereSphereResult},
};

#[test]
fn test_kd_tree() {
    struct LatticeModel {
        lattice_parameters: LatticeVectors,
        atoms: Atoms,
    }
    struct Atoms {
        indices: Vec<usize>,
        symbols: Vec<ElementSymbol>,
        fracs: Vec<CoordData>,
        labels: Vec<Option<String>>,
    }

    impl CoreAtomData for Atoms {
        fn indices(&self) -> &[usize] {
            &self.indices
        }

        fn symbols(&self) -> &[ElementSymbol] {
            &self.symbols
        }

        fn coords(&self) -> &[CoordData] {
            &self.fracs
        }

        fn labels(&self) -> &[Option<String>] {
            &self.labels
        }
    }
    impl CrystalModel for LatticeModel {
        fn get_cell_parameters(&self) -> &impl UnitCellParameters {
            &self.lattice_parameters
        }

        fn get_atom_data(&self) -> &impl CoreAtomData {
            &self.atoms
        }
    }
    let root_dir = env!("CARGO_MANIFEST_DIR");
    let cell_path = Path::new(root_dir).join("SAC_GDY_V.cell");
    let content = read_to_string(cell_path).unwrap();
    let cell_model: CellDocument = CellParser::from_str(&content).parse().unwrap();
    let lattice_param = cell_model.lattice();
    let mut atoms = Atoms {
        indices: Vec::new(),
        symbols: Vec::new(),
        fracs: Vec::new(),
        labels: Vec::new(),
    };
    if let LatticeParam::LatticeCart(lat_cart) = lattice_param {
        let data = [
            lat_cart.a().to_vec(),
            lat_cart.b().to_vec(),
            lat_cart.c().to_vec(),
        ];
        let lattice_vec = LatticeVectors::new(Matrix3::from_vec(data.concat()));
        cell_model
            .ionic_positions()
            .iter()
            .enumerate()
            .for_each(|(id, pos)| {
                atoms.indices.push(id);
                atoms.symbols.push(pos.symbol());
                atoms
                    .fracs
                    .push(CoordData::Fractional(Point3::from(pos.coordinate())));
                atoms.labels.push(None);
            });
        let model = LatticeModel {
            lattice_parameters: lattice_vec,
            atoms,
        };
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
        points.iter().for_each(|p| println!("{}", p));
        let points_tree = build_kd_tree_from_points(&points);
        let dist: f64 = 1.6;
        let mut sphere_sphere_set: HashSet<[usize; 2]> = HashSet::new();
        let mut circle_sphere_set: HashSet<[usize; 3]> = HashSet::new();
        let spheres_results: Vec<Vec<CoordResult>> = points
            .iter()
            .enumerate()
            .map(|(curr, &p)| -> Vec<CoordResult> {
                let query: [f64; 3] = p.into();
                let sphere = Sphere::new(p, dist);
                points_tree
                    .within::<SquaredEuclidean>(&query, 4.0 * dist.powi(2))
                    .iter()
                    .skip(1)
                    .filter_map(|nb| {
                        let nb_id = nb.item as usize;
                        let mut id_pair = [curr, nb_id];
                        id_pair.sort();
                        //
                        if sphere_sphere_set.insert(id_pair) {
                            let nb_sphere = Sphere::new(points[nb_id], dist);
                            match sphere.intersect(&nb_sphere) {
                                SphereSphereResult::Circle(c) => {
                                    let coord_circle = CoordCircle::new(c, [curr, nb_id]);
                                    if let Some(circle_intersect_results) = coord_circle
                                        .nearest_intersect_search(
                                            &points_tree,
                                            &points,
                                            dist,
                                            &mut circle_sphere_set,
                                        )
                                    {
                                        if circle_intersect_results
                                            .iter()
                                            .all(|res| matches!(res, CoordResult::Empty))
                                        {
                                            Some(CoordResult::Circle(coord_circle))
                                        } else {
                                            Some(CoordResult::Various(circle_intersect_results))
                                        }
                                    } else {
                                        Some(CoordResult::Invalid)
                                    }
                                }
                                SphereSphereResult::Point(p) => Some(CoordResult::SinglePoint(
                                    CoordPoint::new(p, vec![curr, nb_id]),
                                )),
                                _ => Some(CoordResult::Empty),
                            }
                        } else {
                            None
                        }
                    })
                    .collect()
            })
            .collect();
        spheres_results.iter().enumerate().for_each(
            |(i, sub_result): (usize, &Vec<CoordResult>)| {
                println!("Atom {i}: ");
                let sphere_sphere_single_points: Vec<&CoordResult> = sub_result
                    .iter()
                    .filter(|res| matches!(res, CoordResult::SinglePoint(p)))
                    .collect();

                sub_result.iter().for_each(|res| match res {
                    CoordResult::Invalid => println!("\tInvalid"),
                    CoordResult::Empty => println!("\tEmpty"),
                    CoordResult::Sphere(_) => println!("\tSphere"),
                    CoordResult::Circle(_) => println!("\tPure Circle"),
                    CoordResult::SinglePoint(p) => println!("\tPoint {p:#?}"),
                    CoordResult::DoublePoints(_, _) => println!("\tDouble point"),
                    CoordResult::Various(results) => {
                        println!("\tVarious:");
                        results.iter().for_each(|r| match r {
                            CoordResult::Invalid => println!("\t\tInvalid"),
                            CoordResult::Empty => (),
                            CoordResult::Sphere(_) => todo!(),
                            CoordResult::Circle(_) => println!("\t\tCircle"),
                            CoordResult::SinglePoint(p) => println!(
                                "\t\tSingle point {}, connects: {:?}",
                                lattice_vec.tensor().try_inverse().unwrap() * p.point,
                                p.atom_ids
                            ),
                            CoordResult::DoublePoints(p1, p2) => {
                                println!(
                                    "\t\tDouble point\n{}, {:?};\n{}, {:?}",
                                    lattice_vec.tensor().try_inverse().unwrap() * p1.point,
                                    p1.atom_ids,
                                    lattice_vec.tensor().try_inverse().unwrap() * p2.point,
                                    p2.atom_ids
                                )
                            }
                            CoordResult::Various(_) => (),
                        })
                    }
                })
            },
        );
    }
}
