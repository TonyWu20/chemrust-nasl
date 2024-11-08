use std::{collections::HashSet, ops::ControlFlow};

use kd_tree::KdIndexTree;
use nalgebra::{distance_squared, Point3};

use crate::geometry::{approx_cmp_f64, approx_eq_point_f64, FloatEq, FloatOrdering};

#[derive(Debug, Clone, PartialEq)]
pub struct MultiCoordPoint {
    pub(crate) point: Point3<f64>,
    pub(crate) atom_ids: Vec<usize>,
}

impl MultiCoordPoint {
    pub fn new(point: Point3<f64>, atom_ids: Vec<usize>) -> Self {
        Self { point, atom_ids }
    }
    pub(crate) fn cn(&self) -> usize {
        self.atom_ids.len()
    }
    pub(crate) fn merge_with(&self, rhs: &Self) -> Option<MultiCoordPoint> {
        if let FloatEq::Eq = approx_eq_point_f64(self.point, rhs.point) {
            let new_connecting_atoms = [self.atom_ids.clone(), rhs.atom_ids.clone()].to_vec();
            let new_connecting_atoms_set: HashSet<usize> =
                new_connecting_atoms.concat().into_iter().collect();
            let mut new_connecting_atoms_array: Vec<usize> =
                new_connecting_atoms_set.into_iter().collect();
            new_connecting_atoms_array.sort();
            Some(MultiCoordPoint::new(self.point, new_connecting_atoms_array))
        } else {
            None
        }
    }

    pub fn no_closer_atoms(
        self,
        kdtree: &KdIndexTree<Point3<f64>>,
        dist: f64,
    ) -> Option<MultiCoordPoint> {
        let closer_than_dist = kdtree
            .within_radius(&self.point(), dist)
            .iter()
            .try_for_each(|&&nb| {
                let distance = distance_squared(&self.point(), kdtree.item(nb));
                if let FloatOrdering::Less = approx_cmp_f64(distance, dist.powi(2)) {
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
        points: &[MultiCoordPoint],
        kdtree: &KdIndexTree<Point3<f64>>,
        dist: f64,
    ) -> Vec<MultiCoordPoint> {
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
        curr: (usize, &MultiCoordPoint),
        points: &[MultiCoordPoint],
        visited: &mut [bool],
    ) -> Option<MultiCoordPoint> {
        let (now, curr_p) = curr;
        if !visited[now] {
            visited[now] = true;
            let same_points: Vec<MultiCoordPoint> = points
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

    pub fn atom_ids(&self) -> &[usize] {
        self.atom_ids.as_ref()
    }

    pub fn point(&self) -> Point3<f64> {
        self.point
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct DelegatePoint<const N: usize> {
    pub(crate) point: Point3<f64>,
    pub(crate) atom_ids: [usize; N],
}

impl<const N: usize> DelegatePoint<N> {
    pub fn new(point: Point3<f64>, atom_ids: [usize; N]) -> Self {
        Self { point, atom_ids }
    }

    pub fn point(&self) -> Point3<f64> {
        self.point
    }

    pub fn atom_ids(&self) -> &[usize] {
        &self.atom_ids
    }
}
