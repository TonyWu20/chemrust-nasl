use std::collections::HashSet;

use kiddo::{ImmutableKdTree, SquaredEuclidean};
use nalgebra::Point3;

use crate::geometry::{
    approx_cmp_f64, approx_eq_point_f64, Circle3d, CircleSphereIntersection, FloatEq,
    FloatOrdering, Intersect, Sphere,
};

#[derive(Debug, Clone, PartialEq)]
pub enum CoordResult {
    Invalid,
    Empty,
    Sphere(CoordSphere),
    Circle(CoordCircle),
    SinglePoint(CoordPoint),
    DoublePoints(CoordPoint, CoordPoint),
    Various(Vec<CoordResult>),
}

impl CoordResult {
    pub fn try_into_sphere(self) -> Result<CoordSphere, Self> {
        if let Self::Sphere(v) = self {
            Ok(v)
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
    pub fn try_into_double_points(self) -> Result<(CoordPoint, CoordPoint), Self> {
        if let Self::DoublePoints(v1, v2) = self {
            Ok((v1, v2))
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
    pub(crate) fn nearest_intersect_search(
        &self,
        kdtree: &ImmutableKdTree<f64, 3>,
        points: &[Point3<f64>],
        dist: f64,
        visited_set: &mut HashSet<[usize; 3]>,
    ) -> Option<Vec<CoordResult>> {
        let neighbors: Vec<Vec<usize>> = self
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
        // Only common neighbors of the associated atoms are possible to
        // form further connections
        let mut common_neighbors: HashSet<usize> = neighbors.concat().iter().cloned().collect();
        self.atom_ids.iter().for_each(|i| {
            common_neighbors.remove(i);
        });
        let neighbor_results: Vec<CoordResult> = common_neighbors
            .iter()
            .filter_map(|&i| {
                let mut id_pairs = [self.atom_ids[0], self.atom_ids[1], i];
                id_pairs.sort();
                if visited_set.insert(id_pairs) {
                    let p = points[i];
                    // circle-sphere intersection
                    let sphere = Sphere::new(p, dist);
                    let circle_sphere = self.circle.intersect(&sphere);
                    match circle_sphere {
                        CircleSphereIntersection::Zero => Some(CoordResult::Empty),
                        CircleSphereIntersection::Single(p) => {
                            Some(CoordResult::SinglePoint(CoordPoint {
                                point: p,
                                atom_ids: [self.atom_ids[0], self.atom_ids[1], i].to_vec(),
                            }))
                        }
                        CircleSphereIntersection::Double(p1, p2) => {
                            let p = if let FloatOrdering::Greater = approx_cmp_f64(p1.z, p2.z) {
                                p1
                            } else {
                                p2
                            };
                            Some(CoordResult::SinglePoint(CoordPoint::new(
                                p,
                                [self.atom_ids[0], self.atom_ids[1], i].to_vec(),
                            )))
                        }
                        // Actually impossible in our current case that every sphere
                        // has the same radius, and thus there can't be a circle has the
                        // same radius as the sphere
                        CircleSphereIntersection::Circle(_) => Some(CoordResult::Invalid),
                        CircleSphereIntersection::InsideSphere => Some(CoordResult::Invalid),
                    }
                } else {
                    None
                }
            })
            .collect();
        if neighbor_results
            .iter()
            .any(|res| matches!(res, CoordResult::Invalid))
        {
            None
        } else {
            Some(neighbor_results)
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
    pub(crate) fn merge_with(self, rhs: Self) -> Option<CoordPoint> {
        if let FloatEq::Eq = approx_eq_point_f64(self.point, rhs.point) {
            let new_connecting_atoms = [self.atom_ids, rhs.atom_ids].to_vec();
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
}
