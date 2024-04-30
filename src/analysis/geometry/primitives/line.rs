use nalgebra::{Point3, UnitVector3};

#[derive(Debug, Clone, Copy)]
pub struct Line {
    origin: Point3<f64>,
    direction: UnitVector3<f64>,
}

impl Line {
    pub fn new(origin: Point3<f64>, direction: UnitVector3<f64>) -> Self {
        Self { origin, direction }
    }

    pub fn origin(&self) -> Point3<f64> {
        self.origin
    }
    pub fn direction(&self) -> UnitVector3<f64> {
        self.direction
    }
    pub fn point(&self, t: f64) -> Point3<f64> {
        self.origin + self.direction().scale(t)
    }
}
