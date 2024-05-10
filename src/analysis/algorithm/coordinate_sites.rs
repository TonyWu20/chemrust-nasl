use nalgebra::Point3;

use crate::analysis::geometry::{Circle3d, Sphere};

#[derive(Debug, Clone, PartialEq)]
pub enum CoordResult {
    Empty,
    Sphere(CoordSphere),
    Circle(CoordCircle),
    Point(CoordPoint),
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
}
