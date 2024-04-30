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
}
