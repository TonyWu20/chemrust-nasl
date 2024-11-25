use std::collections::HashSet;

use nalgebra::Point3;

use crate::geometry::{Intersect, Sphere, SphereSphereResult};

use super::helpers::EnhancedTree;

use crate::coordination_sites::{CoordCircle, CoordResult, MultiCoordPoint};

#[derive(Debug)]
pub struct SphereCheckResult {
    single_points: Vec<MultiCoordPoint>,
    unchecked_circles: Vec<CoordCircle>,
}

impl SphereCheckResult {
    pub fn new(single_points: Vec<MultiCoordPoint>, unchecked_circles: Vec<CoordCircle>) -> Self {
        Self {
            single_points,
            unchecked_circles,
        }
    }

    pub fn single_points(&self) -> &[MultiCoordPoint] {
        self.single_points.as_ref()
    }

    pub fn unchecked_circles(&self) -> &[CoordCircle] {
        self.unchecked_circles.as_ref()
    }
}

pub(crate) fn sphere_check_fn<'a, E>(
    atom_id: usize,
    query: Point3<f64>,
    coord_tree: &'a E,
    radius: f64,
) -> CoordResult
where
    E: EnhancedTree<'a, Result = &'a usize, Item = &'a Point3<f64>>,
{
    let sphere = Sphere::new(query, radius);
    let neighbours = coord_tree.get_sphere_neighbors(query, radius * 2.0);
    if neighbours.len() == 1 {
        CoordResult::Empty
    } else {
        let mut visited_pair: HashSet<[usize; 2]> = HashSet::new();
        let sphere_neighbor_results: Vec<CoordResult> = neighbours
            .iter()
            .skip(1)
            .filter_map(|&&nb_id| {
                let mut id_pair = [atom_id, nb_id];
                id_pair.sort();
                if visited_pair.insert(id_pair) {
                    let nb_sphere = Sphere::new(*coord_tree.item(nb_id), radius);
                    match sphere.intersect(&nb_sphere) {
                        SphereSphereResult::Empty => None,
                        SphereSphereResult::Point(p) => {
                            let coord_point = MultiCoordPoint::new(p, id_pair.to_vec());
                            coord_point
                                .no_closer_atoms(coord_tree, radius)
                                .map(CoordResult::SinglePoint)
                        }
                        SphereSphereResult::Circle(c) => {
                            Some(CoordResult::Circle(CoordCircle::new(c, id_pair)))
                        }
                        SphereSphereResult::Overlap(_) => None,
                    }
                } else {
                    dbg!(id_pair);
                    None
                }
            })
            .collect();
        CoordResult::Various(sphere_neighbor_results)
    }
}
