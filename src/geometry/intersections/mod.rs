mod circle_circle;
mod circle_sphere;
mod plane_plane;
mod sphere_sphere;

use nalgebra::Point3;
use std::f64::EPSILON;

pub use circle_circle::CircleCircleIntersection;
pub use circle_sphere::CircleSphereIntersection;
pub use plane_plane::PlanePlaneIntersection;
pub use sphere_sphere::SphereSphereResult;

#[derive(Debug, Clone, Copy)]
pub enum FloatOrdering {
    Less,
    Equal,
    Greater,
}

#[derive(Debug, Clone, Copy)]
pub enum FloatEq {
    NotEq,
    Eq,
}

pub fn approx_cmp_f64(v1: f64, v2: f64) -> FloatOrdering {
    if v1 - v2 > 5.0 * EPSILON {
        FloatOrdering::Greater
    } else if v1 - v2 < -5.0 * EPSILON {
        FloatOrdering::Less
    } else {
        FloatOrdering::Equal
    }
}

pub fn approx_eq_point_f64(p1: Point3<f64>, p2: Point3<f64>) -> FloatEq {
    let d = p1 - p2;
    // Δx, Δy, Δz < ϵ
    // (Δx^2 + Δy^2 + Δz^2 < 3ϵ^2)
    if d.norm_squared() < 3.0 * (5.0 * EPSILON).powi(2) {
        FloatEq::Eq
    } else {
        FloatEq::NotEq
    }
}

pub trait Intersect<Rhs = Self> {
    type Output;
    fn intersect(&self, rhs: &Rhs) -> Self::Output;
}
