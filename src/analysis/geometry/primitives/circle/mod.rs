use nalgebra::{Point3, UnitVector3};

#[derive(Debug)]
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
}
