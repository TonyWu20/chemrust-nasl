use std::collections::HashSet;

use kiddo::{ImmutableKdTree, SquaredEuclidean};
use nalgebra::Point3;

use crate::{
    geometry::{
        approx_cmp_f64, Circle3d, CircleSphereIntersection, FloatOrdering, Intersect, Sphere,
    },
    CoordPoint, CoordResult,
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
    fn get_common_neighbours(
        &self,
        kdtree: &ImmutableKdTree<f64, 3>,
        points: &[Point3<f64>],
        dist: f64,
    ) -> HashSet<usize> {
        let each_neighbors: Vec<Vec<usize>> = self
            .atom_ids
            .iter()
            .map(|&i| {
                let query: [f64; 3] = points[i].into();
                kdtree
                    .within::<SquaredEuclidean>(&query, 4.0 * dist.powi(2))
                    .iter()
                    .skip(1)
                    .map(|nb| nb.item as usize)
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
            let coord_points: Vec<CoordPoint> = neighbour_results
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
        kdtree: &ImmutableKdTree<f64, 3>,
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

                // #[cfg(debug_assertions)]
                // {
                //     if i == 44 && self.atom_ids() == [24, 26] {
                //         let cs_cc = self.circle().center() - sphere.center();
                //         let cut_at = self.circle().n().dot(&cs_cc);
                //         println!("Cut at and radius");
                //         dbg!(cut_at.abs(), sphere.radius());
                //         let new_circle_center = p + self.circle().n().scale(cut_at);
                //         let new_circle_radius = (dist.powi(2) - cut_at.powi(2)).sqrt();
                //         dbg!(new_circle_center, new_circle_radius);
                //         dbg!(self.circle.center(), self.circle.radius());
                //         let c1c2 = new_circle_center - self.circle.center();
                //         let r1r2_sum_squared = (new_circle_radius + self.circle.radius()).powi(2);
                //         dbg!(r1r2_sum_squared - c1c2.norm_squared());
                //         dbg!(r1r2_sum_squared.sqrt() - c1c2.norm());
                //         dbg!(approx_cmp_f64(c1c2.norm_squared(), r1r2_sum_squared));
                //         let r1r2_diff_squared = (new_circle_radius - self.circle.radius()).powi(2);
                //         dbg!(approx_cmp_f64(c1c2.norm_squared(), r1r2_diff_squared));
                //         dbg!(circle_sphere);
                //     }
                // }
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
                CoordResult::SinglePoint(CoordPoint::new(p, atom_id))
            }
            CircleSphereIntersection::Double(p1, p2) => {
                let p = if let FloatOrdering::Greater = approx_cmp_f64(p1.z, p2.z) {
                    p1
                } else {
                    p2
                };
                let mut atom_id = [circle_id[0], circle_id[1], sphere_id].to_vec();
                atom_id.sort();
                CoordResult::SinglePoint(CoordPoint::new(p, atom_id))
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
