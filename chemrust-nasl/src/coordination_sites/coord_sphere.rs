use crate::geometry::Sphere;

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct CoordSphere {
    pub(crate) sphere: Sphere,
    pub(crate) atom_id: usize,
}

impl CoordSphere {
    pub fn new(sphere: Sphere, atom_id: usize) -> Self {
        Self { sphere, atom_id }
    }

    pub fn atom_id(&self) -> usize {
        self.atom_id
    }

    pub fn sphere(&self) -> Sphere {
        self.sphere
    }
}
