use nalgebra::{Point3, UnitVector3, Vector3};

#[derive(Debug, Clone, Copy)]
/// n1*x + n2*y + n3*z + d = 0
/// # Note
/// The normal is a unit vector
/// So if two planes have different d and cross product of two norms is a 0 vector,
/// they must be parallel to each other.
pub struct Plane {
    normal: UnitVector3<f64>,
    d: f64,
}

impl Plane {
    pub fn new(normal: UnitVector3<f64>, d: f64) -> Self {
        Self { normal, d }
    }
    pub fn from_normal_and_point(normal: UnitVector3<f64>, point: Point3<f64>) -> Self {
        let op: Vector3<f64> = point - Point3::origin();
        let d = normal.dot(&op) * -1.0;
        Self { normal, d }
    }

    pub fn normal(&self) -> UnitVector3<f64> {
        self.normal
    }
    pub fn d(&self) -> f64 {
        self.d
    }
}
