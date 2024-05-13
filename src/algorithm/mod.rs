use std::f64::EPSILON;

use kiddo::{ImmutableKdTree, NearestNeighbour, SquaredEuclidean};
use nalgebra::Point3;

use super::algorithm::coordinate_sites::CoordSphere;

use self::coordinate_sites::{CoordCircle, CoordPoint, CoordResult};

use super::geometry::{Intersect, Sphere};

mod circle_check;
mod coordinate_sites;
mod sphere_check;

#[derive(Debug, Clone, Copy)]
pub struct Config<'a> {
    to_check: &'a [Point3<f64>],
    bondlength: f64,
}

#[derive(Debug)]
pub struct SiteKdTree(ImmutableKdTree<f64, 3>);

#[cfg(test)]
mod test;

fn build_kd_tree_from_points(points: &[Point3<f64>]) -> ImmutableKdTree<f64, 3> {
    let entries = points.iter().map(|&p| p.into()).collect::<Vec<[f64; 3]>>();
    ImmutableKdTree::new_from_slice(&entries)
}

fn build_kd_tree_from_circles(circles: &[CoordCircle]) -> ImmutableKdTree<f64, 3> {
    let entries: Vec<[f64; 3]> = circles
        .iter()
        .map(|&circ| circ.circle.center().into())
        .collect();
    ImmutableKdTree::new_from_slice(&entries)
}

pub struct SiteFinder {
    coords: Vec<Point3<f64>>,
    coord_tree: ImmutableKdTree<f64, 3>,
    bondlength: f64,
}

impl SiteFinder {
    pub fn new_from_points(points: &[Point3<f64>], bondlength: f64) -> Self {
        SiteFinder {
            coords: points.to_vec(),
            coord_tree: build_kd_tree_from_points(points),
            bondlength,
        }
    }
    fn build_sphere_from_id(&self, id: usize) -> Sphere {
        Sphere::new(self.coords[id], self.bondlength)
    }
    /// Search with the given bondlength.
    /// # Steps
    /// 0. A filter by xyz range could be possibly conducted
    /// 1. Iterate every selected coord, search the nearest neighbor within the radius of 2 x bondlength (minimum requirement for sphere-sphere intersection) by the KdTree constructed with all coordinates
    /// 2. For each result of NN:
    /// - Only the coord itself: The new site is on the surface of a sphere with a center on this coord, and a radius of the bondlength.
    /// - More than itself:
    /// - check sphere-sphere intersect of coord itself with every other coord, possible outcomes:
    /// - A `CoordPoint`
    /// - A `Circle3d`
    /// - Have a `Vec<SphereSphereResult>` for each result
    /// - Convert into `Vec<CoordResult>`
    /// 3. Collect all sphere-sphere results, concat into one `Vec<CoordResult>`
    /// 4. Filter `CoordCircle`, check intersect (O(n^2)), collect into another `Vec<CoordResult>`
    /// 5. Filter `CoordPoint` from step **3** and **4**, check duplicate, merge with collecting all relevant atom id
    /// 6. Now we have results of `CoordSphere` from **3**, `CoordCircle` from **4**
    /// ## Recheck
    /// - Circle recheck:
    /// 1. Get nearest neighbors of every circle center within a radius of `circle_r + bondlength`, exclude the two neighbors connected by the circle
    /// 2. Determine the longest possible distance from each atom to the circle
    /// - Compute the distance of atom coord to the plane of the circle
    /// - Connect the coord to center of radius, compute `x = sqrt(l^2 - h^2)`
    /// - Longest possible distance `L^2 = h^2 + x^2 < Bondlength^2`
    /// 3. If `L^2 < Bondlength^2`, remove this circle
    /// 4. else, this circle is out of reach of other points
    /// - Point recheck
    /// 1. Get nearest neighbors of every points within a radius of bondlength + EPSILON
    /// 2. If there is one atom or more that has distance < bondlength, drop
    pub fn search(&self) {
        todo!()
        // let sphere_intersect_results: Vec<Vec<CoordResult>> = self
        //     .coord_tree
        //     .iter()
        //     .map(|(id, coord)| self.sphere_stage(id as usize, &coord))
        //     .collect();
        // let all_coord_results = sphere_intersect_results.concat();
        // let circles: Vec<CoordCircle> = all_coord_results
        //     .iter()
        //     .filter_map(|item| {
        //         if let CoordResult::Circle(coord_circ) = &&item {
        //             Some(*coord_circ)
        //         } else {
        //             None
        //         }
        //     })
        //     .collect();
    }
    // fn sphere_stage(&self, id: usize, coord: &[f64; 3]) -> Vec<CoordResult> {
    //     let sphere_self = self.build_sphere_from_id(id);
    //     let dist = (2.0 * (self.bondlength + EPSILON)).powi(2);
    //     let neighbours = self.coord_tree.within::<SquaredEuclidean>(coord, dist);
    //     if neighbours.len() == 1 {
    //         vec![CoordResult::Sphere(CoordSphere::new(sphere_self, id))]
    //     } else {
    //         self.neighbour_spheres_intersect(&sphere_self, id, &neighbours)
    //     }
    // }
    // fn neighbour_spheres_intersect(
    //     &self,
    //     sphere_self: &Sphere,
    //     id: usize,
    //     neighbours: &[NearestNeighbour<f64, u64>],
    // ) -> Vec<CoordResult> {
    //     neighbours
    //         .iter()
    //         .skip(1)
    //         .map(|neighbour| {
    //             let nb_id = neighbour.item as usize;
    //             let nb_sphere = self.build_sphere_from_id(nb_id);
    //             let result = sphere_self.intersect(&nb_sphere);
    //             match result {
    //                 super::geometry::SphereSphereResult::Empty => CoordResult::Empty,
    //                 super::geometry::SphereSphereResult::Point(p) => {
    //                     let coord_pt = CoordPoint::new(p, vec![id, nb_id]);
    //                     CoordResult::Point(coord_pt)
    //                 }
    //                 super::geometry::SphereSphereResult::Circle(c) => {
    //                     let coord_circ = CoordCircle::new(c, [id, nb_id]);
    //                     CoordResult::Circle(coord_circ)
    //                 }
    //                 super::geometry::SphereSphereResult::Overlap(_) => CoordResult::Empty,
    //             }
    //         })
    //         .collect()
    // }
}
