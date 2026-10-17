use std::cmp::Ordering;

use kd_tree::KdIndexTree;
use nalgebra::{distance_squared, Point3};

use crate::approx_cmp_f64;

pub trait EnhancedTree<'a> {
    type Result;
    type Item;
    fn within_radius_dist_sorted(&'a self, query: Point3<f64>, dist: f64) -> Vec<Self::Result>;
    fn get_sphere_neighbors(&'a self, query: Point3<f64>, dist: f64) -> Vec<Self::Result>;
    fn item(&'a self, id: usize) -> Self::Item;
}

impl<'a> EnhancedTree<'a> for KdIndexTree<'a, Point3<f64>> {
    type Result = &'a usize;
    type Item = &'a Point3<f64>;

    fn within_radius_dist_sorted(&'a self, query: Point3<f64>, dist: f64) -> Vec<Self::Result> {
        let mut neighbours = self.within_radius(&query, dist);
        neighbours.sort_by(|&&a, &&b| {
            let a_to_q = distance_squared(&query, self.item(a));
            let b_to_q = distance_squared(&query, self.item(b));
            match approx_cmp_f64(a_to_q, b_to_q) {
                crate::FloatOrdering::Less => Ordering::Less,
                crate::FloatOrdering::Equal => Ordering::Equal,
                crate::FloatOrdering::Greater => Ordering::Greater,
            }
        });
        neighbours
    }

    fn get_sphere_neighbors(&'a self, query: Point3<f64>, dist: f64) -> Vec<Self::Result> {
        self.within_radius_dist_sorted(query, dist * 2.0)
    }

    fn item(&'a self, id: usize) -> Self::Item {
        self.item(id)
    }
}
