use std::collections::HashSet;

use kiddo::SquaredEuclidean;

use crate::geometry::{Intersect, Sphere, SphereSphereResult};

use super::{SearchConfig, SiteIndex};

use crate::coordination_sites::{CoordCircle, CoordResult, MultiCoordPoint};

#[derive(Debug)]
pub(crate) struct SphereCheckResult {
    single_points: Vec<MultiCoordPoint>,
    unchecked_circles: Vec<CoordCircle>,
    // spheres: Vec<CoordSphere>,
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

    // pub fn spheres(&self) -> &[CoordSphere] {
    //     self.spheres.as_ref()
    // }
}

/// The first round search, abstract into sphere-sphere intersection.
/// If the sphere does not have possible intersecting neighbours, then
/// return early as `CoordResult::Sphere`. Otherwise, Use `CoordResult::Various`
/// to unify the possible `CoordPoint` and `CoordCircle` (cut and intersect of two spheres)
pub fn sphere_check(site_index: &SiteIndex, search_config: &SearchConfig) -> SphereCheckResult {
    let to_check = search_config.to_check;
    let dist = search_config.bondlength;
    let mut visited_pair: HashSet<[usize; 2]> = HashSet::new();
    let mut results: Vec<CoordResult> = to_check
        .iter()
        .map(
            // Use `CoordResult::Various` to unify points and circles
            |&(atom_id, p)| -> CoordResult {
                let query: [f64; 3] = p.into();
                let sphere = Sphere::new(p, dist);
                let sphere_neighbours = site_index
                    .coord_tree
                    .within::<SquaredEuclidean>(&query, 4.0 * dist.powi(2));
                if sphere_neighbours.len() == 1 {
                    CoordResult::Empty
                } else {
                    let sphere_neighbor_results: Vec<CoordResult> = sphere_neighbours
                        .iter()
                        .skip(1)
                        .filter_map(|nb| {
                            let nb_id = nb.item as usize;
                            let mut id_pair = [atom_id, nb_id];
                            id_pair.sort();
                            if visited_pair.insert(id_pair) {
                                let nb_sphere = Sphere::new(site_index.coords()[nb_id], dist);
                                match sphere.intersect(&nb_sphere) {
                                    SphereSphereResult::Empty => None,
                                    SphereSphereResult::Point(p) => {
                                        let coord_point = MultiCoordPoint::new(p, id_pair.to_vec());
                                        coord_point
                                            .no_closer_atoms(
                                                site_index.coord_tree(),
                                                search_config.bondlength,
                                            )
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
            },
        )
        .collect();
    // let spheres: Vec<CoordSphere> = results
    //     .iter()
    //     .filter_map(|res| res.try_into_sphere().ok())
    //     .collect();
    let unchecked_circles: Vec<Vec<CoordCircle>> = results
        .iter()
        .filter_map(|res| res.try_pull_circles_from_various().ok())
        .collect();
    let points: Vec<Vec<MultiCoordPoint>> = results
        .iter_mut()
        .filter_map(|res| res.try_pull_single_points_from_various().ok())
        .collect();
    // SphereCheckResult::new(points.concat(), unchecked_circles.concat(), spheres)
    SphereCheckResult::new(points.concat(), unchecked_circles.concat())
}
