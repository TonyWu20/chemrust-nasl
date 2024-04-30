mod circle_circle;
mod plane_plane;
mod sphere_sphere;

use std::f64::EPSILON;

pub use sphere_sphere::{sphere_sphere_intersect, SphereSphereResult};

pub enum FloatOrdering {
    Less,
    Equal,
    Greater,
}

fn approx_cmp_f64(v1: f64, v2: f64) -> FloatOrdering {
    if v1 - v2 > EPSILON {
        FloatOrdering::Greater
    } else if v1 - v2 < -1.0 * EPSILON {
        FloatOrdering::Less
    } else {
        FloatOrdering::Equal
    }
}
