use nalgebra::{Point3, UnitVector3, Vector3};

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

    pub fn get_point_on_circle(&self, theta: f64) -> Point3<f64> {
        let (x, y, _z) = (self.n().x, self.n().y, self.n().z);
        // We want the v2 to act as the "z-axis" after transformation
        let pre_v1 = Vector3::new(-1.0 * y, x, 0.0);
        let v1 = if pre_v1.norm_squared() < f64::EPSILON {
            // such v1 becomes a null vector when the normal is (0, 0, z)
            Vector3::x_axis()
        } else {
            UnitVector3::new_normalize(pre_v1)
        };
        let v2 = self.n().cross(&v1);
        self.center() + (v1.scale(theta.cos()) + v2.scale(theta.sin())).scale(self.radius())
    }
}
