use std::collections::HashSet;

use kd_tree::KdIndexTree;
use nalgebra::Point3;

use crate::geometry::{Intersect, Sphere, SphereSphereResult};

use super::{SearchConfig, SiteIndex};

use crate::coordination_sites::{CoordCircle, CoordResult, MultiCoordPoint};

#[derive(Debug)]
pub(crate) struct SphereCheckResult {
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

/// The first round search, abstract into sphere-sphere intersection.
/// If the sphere does not have possible intersecting neighbours, then
/// return early as `CoordResult::Sphere`. Otherwise, Use `CoordResult::Various`
/// to unify the possible `CoordPoint` and `CoordCircle` (cut and intersect of two spheres)
pub fn sphere_check(site_index: &SiteIndex, search_config: &SearchConfig) -> SphereCheckResult {
    let to_check = search_config.to_check;
    let mut results: Vec<CoordResult> = to_check
        .iter()
        .map(
            // Use `CoordResult::Various` to unify points and circles
            |&(atom_id, p)| -> CoordResult {
                sphere_check_fn(atom_id, p, site_index, search_config)
            },
        )
        .collect();
    let unchecked_circles: Vec<Vec<CoordCircle>> = results
        .iter()
        .filter_map(|res| res.try_pull_circles_from_various().ok())
        .collect();
    let points: Vec<Vec<MultiCoordPoint>> = results
        .iter_mut()
        .filter_map(|res| res.try_pull_single_points_from_various().ok())
        .collect();
    SphereCheckResult::new(points.concat(), unchecked_circles.concat())
}

fn sphere_check_fn(
    atom_id: usize,
    query: Point3<f64>,
    site_index: &SiteIndex,
    search_config: &SearchConfig,
) -> CoordResult {
    let sphere = Sphere::new(query, search_config.bondlength());
    let tree = site_index.coord_tree();
    let neighbours = tree.within_radius(&query, search_config.bondlength());
    if neighbours.len() == 1 {
        CoordResult::Empty
    } else {
        sphere_neighbour_check(
            &sphere,
            atom_id,
            &neighbours,
            site_index.coord_tree(),
            search_config.bondlength(),
        )
    }
}

fn sphere_neighbour_check(
    sphere: &Sphere,
    atom_id: usize,
    neighbours: &[&usize],
    coord_tree: &KdIndexTree<Point3<f64>>,
    dist: f64,
) -> CoordResult {
    let mut visited_pair: HashSet<[usize; 2]> = HashSet::new();
    let sphere_neighbor_results: Vec<CoordResult> = neighbours
        .iter()
        .skip(1)
        .filter_map(|&&nb_id| {
            let mut id_pair = [atom_id, nb_id];
            id_pair.sort();
            if visited_pair.insert(id_pair) {
                let nb_sphere = Sphere::new(*coord_tree.item(nb_id), dist);
                match sphere.intersect(&nb_sphere) {
                    SphereSphereResult::Empty => None,
                    SphereSphereResult::Point(p) => {
                        let coord_point = MultiCoordPoint::new(p, id_pair.to_vec());
                        coord_point
                            .no_closer_atoms(coord_tree, dist)
                            .map(CoordResult::SinglePoint)
                    }
                    SphereSphereResult::Circle(c) => {
                        Some(CoordResult::Circle(CoordCircle::new(c, id_pair)))
                    }
                    SphereSphereResult::Overlap(_) => None,
                }
            } else {
                None
            }
        })
        .collect();
    CoordResult::Various(sphere_neighbor_results)
}
