use nalgebra::{Point3, UnitVector3};

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Sphere {
    center: Point3<f64>,
    radius: f64,
}

impl Sphere {
    pub fn new(center: Point3<f64>, radius: f64) -> Self {
        Self { center, radius }
    }
    pub fn point_at_surface(&self, direction: &UnitVector3<f64>) -> Point3<f64> {
        self.center + direction.scale(self.radius)
    }

    pub fn center(&self) -> Point3<f64> {
        self.center
    }

    pub fn radius(&self) -> f64 {
        self.radius
    }
}
