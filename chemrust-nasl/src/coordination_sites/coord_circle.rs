use std::{
    collections::HashSet,
    f64::consts::{FRAC_PI_2, FRAC_PI_8},
    ops::ControlFlow,
};

use kd_tree::KdIndexTree;
use nalgebra::Point3;

use crate::{
    algorithm::EnhancedTree,
    geometry::{
        approx_cmp_f64, Circle3d, CircleSphereIntersection, FloatOrdering, Intersect, Sphere,
    },
    CoordResult, DelegatePoint, MultiCoordPoint,
};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CoordCircle {
    pub(crate) circle: Circle3d,
    pub(crate) atom_ids: [usize; 2],
}

impl CoordCircle {
    pub fn new(circle: Circle3d, atom_ids: [usize; 2]) -> Self {
        Self { circle, atom_ids }
    }

    pub fn get_possible_point(
        &self,
        coord_tree: &KdIndexTree<Point3<f64>>,
        dist: f64,
    ) -> Option<DelegatePoint<2>> {
        fn each_possible_theta(step: i32) -> f64 {
            let step_frac_pi_32 = FRAC_PI_8 / 4.0;
            FRAC_PI_2 + step as f64 * step_frac_pi_32
        }
        let possible_position = (0..32).map(each_possible_theta).try_for_each(|theta| {
            let query = self.circle().get_point_on_circle(theta);
            let neighbours = coord_tree.within_radius(&query, dist + 1e-5_f64);
            if neighbours.len() == 2
                && neighbours.contains(&&self.atom_ids()[0])
                && neighbours.contains(&&self.atom_ids()[1])
            {
                ControlFlow::Break(query)
            } else {
                ControlFlow::Continue(())
            }
        });
        match possible_position {
            ControlFlow::Continue(_) => None,
            ControlFlow::Break(pos) => Some(DelegatePoint::<2>::new(pos, self.atom_ids)),
        }
    }

    fn get_common_neighbours(
        &self,
        kdtree: &KdIndexTree<Point3<f64>>,
        points: &[Point3<f64>],
        dist: f64,
    ) -> HashSet<usize> {
        let each_neighbors: Vec<Vec<usize>> = self
            .atom_ids
            .iter()
            .map(|&i| {
                let query = points[i];
                kdtree
                    .within_radius_dist_sorted(query, 2.0 * (dist + 1e-5_f64))
                    .iter()
                    .skip(1)
                    .map(|&&i| i)
                    .collect()
            })
            .collect();
        let mut common_neighbors: HashSet<usize> =
            each_neighbors.concat().iter().cloned().collect();
        // Remove self ids included when searching NNs.
        self.atom_ids.iter().for_each(|i| {
            common_neighbors.remove(i);
        });
        common_neighbors
    }
    fn classify_neighbour_results(
        &self,
        neighbour_results: Vec<CoordResult>,
    ) -> Option<CoordResult> {
        if neighbour_results
            .iter()
            .any(|res| matches!(res, CoordResult::Invalid))
        {
            None
        } else if neighbour_results
            .iter()
            .all(|res| matches!(res, CoordResult::Empty))
        {
            Some(CoordResult::Circle(*self))
        } else {
            let coord_points: Vec<MultiCoordPoint> = neighbour_results
                .iter()
                .filter_map(|res| {
                    if let CoordResult::SinglePoint(coord_point) = res {
                        Some(coord_point.clone())
                    } else {
                        None
                    }
                })
                .collect();
            Some(CoordResult::Points(coord_points))
        }
    }
    /// # Returns
    /// `None` when there is atoms closer than the required distance
    /// `Some` for 1. `CoordResult::Circle`2. `CoordResult::Points`
    pub(crate) fn common_neighbours_intersect(
        &self,
        kdtree: &KdIndexTree<Point3<f64>>,
        points: &[Point3<f64>],
        dist: f64,
    ) -> Option<CoordResult> {
        // Only common neighbors of the associated atoms are possible to
        // form further connections
        let common_neighbors: HashSet<usize> = self.get_common_neighbours(kdtree, points, dist);
        let neighbor_results: Vec<CoordResult> = common_neighbors
            .iter()
            .map(|&i| {
                let p = points[i];
                // circle-sphere intersection
                let sphere = Sphere::new(p, dist);
                let circle_sphere = self.circle.intersect(&sphere);
                circle_sphere.to_coord_result(&self.atom_ids, i)
            })
            .collect();
        self.classify_neighbour_results(neighbor_results)
    }

    pub fn atom_ids(&self) -> [usize; 2] {
        self.atom_ids
    }

    pub fn circle(&self) -> Circle3d {
        self.circle
    }
}

impl CircleSphereIntersection {
    pub fn to_coord_result(self, circle_id: &[usize; 2], sphere_id: usize) -> CoordResult {
        match self {
            CircleSphereIntersection::Zero => CoordResult::Empty,
            CircleSphereIntersection::Single(p) => {
                let mut atom_id = [circle_id[0], circle_id[1], sphere_id].to_vec();
                atom_id.sort();
                CoordResult::SinglePoint(MultiCoordPoint::new(p, atom_id))
            }
            CircleSphereIntersection::Double(p1, p2) => {
                let p = if let FloatOrdering::Greater = approx_cmp_f64(p1.z, p2.z) {
                    p1
                } else {
                    p2
                };
                let mut atom_id = [circle_id[0], circle_id[1], sphere_id].to_vec();
                atom_id.sort();
                CoordResult::SinglePoint(MultiCoordPoint::new(p, atom_id))
            }
            // Actually impossible in our current case that every sphere
            // has the same radius, and thus there can't be a circle has the
            // same radius as the sphere
            CircleSphereIntersection::Circle(_) => CoordResult::Invalid,
            CircleSphereIntersection::InsideSphere => CoordResult::Invalid,
            CircleSphereIntersection::SphereInCircle => CoordResult::Invalid,
            CircleSphereIntersection::Invalid => CoordResult::Invalid,
        }
    }
}
