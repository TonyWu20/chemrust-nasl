use std::{collections::HashSet, ops::ControlFlow};

use kiddo::{ImmutableKdTree, SquaredEuclidean};
use nalgebra::Point3;

use crate::geometry::{
    approx_cmp_f64, approx_eq_point_f64, Circle3d, CircleSphereIntersection, FloatEq,
    FloatOrdering, Intersect, Sphere,
};

mod visualize;

pub use visualize::*;

#[derive(Debug, Clone, PartialEq)]
pub enum CoordResult {
    Invalid,
    Empty,
    Sphere(CoordSphere),
    Circle(CoordCircle),
    SinglePoint(CoordPoint),
    Points(Vec<CoordPoint>),
    Various(Vec<CoordResult>),
}

impl CoordResult {
    pub fn try_into_sphere(&self) -> Result<CoordSphere, &Self> {
        if let Self::Sphere(v) = self {
            Ok(*v)
        } else {
            Err(self)
        }
    }

    pub fn try_into_circle(self) -> Result<CoordCircle, Self> {
        if let Self::Circle(v) = self {
            Ok(v)
        } else {
            Err(self)
        }
    }

    pub fn try_into_single_point(self) -> Result<CoordPoint, Self> {
        if let Self::SinglePoint(v) = self {
            Ok(v)
        } else {
            Err(self)
        }
    }

    pub fn try_into_points(self) -> Result<Vec<CoordPoint>, Self> {
        if let Self::Points(v) = self {
            Ok(v)
        } else {
            Err(self)
        }
    }
    pub fn try_pull_single_points_from_various(&self) -> Result<Vec<CoordPoint>, &Self> {
        if let Self::Various(v) = self {
            Ok(v.iter()
                .filter_map(|res| {
                    if let CoordResult::SinglePoint(p) = res {
                        Some(p.to_owned())
                    } else {
                        None
                    }
                })
                .collect())
        } else {
            Err(self)
        }
    }
    pub fn try_pull_points_from_various(&self) -> Result<Vec<CoordPoint>, &Self> {
        if let Self::Various(v) = self {
            let points: Vec<Vec<CoordPoint>> = v
                .iter()
                .filter_map(|res| {
                    if let CoordResult::Points(p) = res {
                        Some(p.to_owned())
                    } else {
                        None
                    }
                })
                .collect();
            Ok(points.concat())
        } else {
            Err(self)
        }
    }
    pub fn try_pull_circles_from_various(&self) -> Result<Vec<CoordCircle>, &Self> {
        if let Self::Various(v) = self {
            let circles: Vec<CoordCircle> = v
                .iter()
                .filter_map(|res| {
                    if let CoordResult::Circle(c) = res {
                        Some(*c)
                    } else {
                        None
                    }
                })
                .collect();
            Ok(circles)
        } else {
            Err(self)
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub(crate) struct CoordSphere {
    pub(crate) sphere: Sphere,
    pub(crate) atom_id: usize,
}

impl CoordSphere {
    pub(crate) fn new(sphere: Sphere, atom_id: usize) -> Self {
        Self { sphere, atom_id }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct CoordCircle {
    pub(crate) circle: Circle3d,
    pub(crate) atom_ids: [usize; 2],
}

impl CoordCircle {
    pub(crate) fn new(circle: Circle3d, atom_ids: [usize; 2]) -> Self {
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
                circle_sphere.to_coord_result(&self.atom_ids, i)
            })
            .collect();
        self.classify_neighbour_results(neighbor_results)
    }
}

impl CircleSphereIntersection {
    pub fn to_coord_result(self, circle_id: &[usize; 2], sphere_id: usize) -> CoordResult {
        match self {
            CircleSphereIntersection::Zero => CoordResult::Empty,
            CircleSphereIntersection::Single(p) => CoordResult::SinglePoint(CoordPoint::new(
                p,
                [circle_id[0], circle_id[1], sphere_id].to_vec(),
            )),
            CircleSphereIntersection::Double(p1, p2) => {
                let p = if let FloatOrdering::Greater = approx_cmp_f64(p1.z, p2.z) {
                    p1
                } else {
                    p2
                };
                CoordResult::SinglePoint(CoordPoint::new(
                    p,
                    [circle_id[0], circle_id[1], sphere_id].to_vec(),
                ))
            }
            // Actually impossible in our current case that every sphere
            // has the same radius, and thus there can't be a circle has the
            // same radius as the sphere
            CircleSphereIntersection::Circle(_) => CoordResult::Invalid,
            CircleSphereIntersection::InsideSphere => CoordResult::Invalid,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CoordPoint {
    pub(crate) point: Point3<f64>,
    pub(crate) atom_ids: Vec<usize>,
}

impl CoordPoint {
    pub(crate) fn new(point: Point3<f64>, atom_ids: Vec<usize>) -> Self {
        Self { point, atom_ids }
    }
    pub(crate) fn cn(&self) -> usize {
        self.atom_ids.len()
    }
    pub(crate) fn merge_with(&self, rhs: &Self) -> Option<CoordPoint> {
        if let FloatEq::Eq = approx_eq_point_f64(self.point, rhs.point) {
            let new_connecting_atoms = [self.atom_ids.clone(), rhs.atom_ids.clone()].to_vec();
            let new_connecting_atoms_set: HashSet<usize> =
                new_connecting_atoms.concat().into_iter().collect();
            let mut new_connecting_atoms_array: Vec<usize> =
                new_connecting_atoms_set.into_iter().collect();
            new_connecting_atoms_array.sort();
            Some(CoordPoint::new(self.point, new_connecting_atoms_array))
        } else {
            None
        }
    }

    pub fn no_closer_atoms(
        self,
        kdtree: &ImmutableKdTree<f64, 3>,
        dist: f64,
    ) -> Option<CoordPoint> {
        let query: [f64; 3] = self.point.into();
        let closer_than_dist = kdtree
            .within::<SquaredEuclidean>(&query, dist.powi(2))
            .into_iter()
            .try_for_each(|nb| {
                if let FloatOrdering::Less = approx_cmp_f64(nb.distance, dist.powi(2)) {
                    ControlFlow::Break(nb)
                } else {
                    ControlFlow::Continue(())
                }
            });
        if !matches!(closer_than_dist, ControlFlow::Break(_)) {
            Some(self)
        } else {
            None
        }
    }
    pub fn dedup_points(
        points: &[CoordPoint],
        kdtree: &ImmutableKdTree<f64, 3>,
        dist: f64,
    ) -> Vec<CoordPoint> {
        let mut visited = vec![false; points.len()];
        points
            .iter()
            .enumerate()
            .filter_map(|(now, curr_p)| {
                Self::look_for_same_points((now, curr_p), points, &mut visited)
            })
            .filter_map(|p| p.no_closer_atoms(kdtree, dist))
            .collect()
    }
    fn look_for_same_points(
        curr: (usize, &CoordPoint),
        points: &[CoordPoint],
        visited: &mut [bool],
    ) -> Option<CoordPoint> {
        let (now, curr_p) = curr;
        if !visited[now] {
            visited[now] = true;
            let same_points: Vec<CoordPoint> = points
                .iter()
                .enumerate()
                .filter_map(|(to_check, p)| {
                    if visited[to_check] {
                        None
                    } else if let Some(merged) = curr_p.merge_with(p) {
                        visited[to_check] = true;
                        Some(merged)
                    } else {
                        None
                    }
                })
                .collect();
            if same_points.is_empty() {
                Some(curr_p.clone())
            } else {
                Some(
                    same_points
                        .into_iter()
                        .reduce(|acc, x| acc.merge_with(&x).unwrap())
                        .unwrap(),
                )
            }
        } else {
            None
        }
    }
}
