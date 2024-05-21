use nalgebra::Point3;

use crate::geometry::Sphere;

use super::{approx_cmp_f64, Intersect};

#[derive(Debug, Clone, Copy)]
pub enum PointToSphere {
    Inside,
    Outside,
    OnSurface,
}

impl Intersect<Point3<f64>> for Sphere {
    type Output = PointToSphere;

    fn intersect(&self, rhs: &Point3<f64>) -> Self::Output {
        let op = rhs - self.center();
        match approx_cmp_f64(op.norm_squared(), self.radius().powi(2)) {
            super::FloatOrdering::Less => PointToSphere::Inside,
            super::FloatOrdering::Equal => PointToSphere::OnSurface,
            super::FloatOrdering::Greater => PointToSphere::Outside,
        }
    }
}
