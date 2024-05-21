use nalgebra::{Point3, UnitVector3};

use super::plane::Plane;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Circle3d {
    center: Point3<f64>,
    radius: f64,
    /// Plane normal vector
    n: UnitVector3<f64>,
}

impl Circle3d {
    pub fn new(center: Point3<f64>, radius: f64, n: UnitVector3<f64>) -> Self {
        Self { center, radius, n }
    }

    pub fn center(&self) -> Point3<f64> {
        self.center
    }

    pub fn radius(&self) -> f64 {
        self.radius
    }

    pub fn n(&self) -> UnitVector3<f64> {
        self.n
    }
    pub fn plane_of_circle(&self) -> Plane {
        Plane::from_normal_and_point(self.n(), self.center())
    }
    pub fn point_to_circle_distance_range(&self, point: Point3<f64>) -> (f64, f64) {
        let op = point - self.center;
        let op_normal_angle = self.n.angle(&op);
        let op_projection_distance = op.norm() * op_normal_angle.sin();
        let op_projection_height = op.norm() * op_normal_angle.cos();
        let min_x = (op_projection_distance - self.radius).abs();
        let max_x = op_projection_distance + self.radius;
        let min_dist = (min_x.powi(2) + op_projection_height.powi(2)).sqrt();
        let max_dist = (max_x.powi(2) + op_projection_height.powi(2)).sqrt();
        (min_dist, max_dist)
    }
}
