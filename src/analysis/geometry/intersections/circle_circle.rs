use std::f64::EPSILON;

use nalgebra::{Point3, UnitVector3};

use crate::analysis::geometry::{
    intersections::{approx_cmp_f64, FloatOrdering},
    primitives::Circle3d,
};

#[derive(Debug, Clone, Copy)]
enum CircleCircleRelationship {
    Coplanar,
    ParallelPlane,
    NotParallel,
}

impl CircleCircleRelationship {
    fn determine(c1: &Circle3d, c2: &Circle3d) -> Self {
        let n3 = c1.n().cross(&c2.n());
        if n3.norm_squared() < EPSILON {
            let d = c1.center() - c2.center();
            let d_dot_n1 = d.dot(&c1.n());
            if d_dot_n1.abs() > EPSILON {
                Self::ParallelPlane
            } else {
                Self::Coplanar
            }
        } else {
            Self::NotParallel
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub enum CircleCircleIntersection {
    Empty,
    Single(Point3<f64>),
    Double(Point3<f64>, Point3<f64>),
}

impl CircleCircleIntersection {
    pub fn intersect_result(c1: &Circle3d, c2: &Circle3d) -> Self {
        let case = CircleCircleRelationship::determine(c1, c2);
        match case {
            CircleCircleRelationship::ParallelPlane => CircleCircleIntersection::Empty,
            CircleCircleRelationship::Coplanar => coplanar_circle_circle_intersect(c1, c2),
            CircleCircleRelationship::NotParallel => noncoplanar_circle_circle_intersect(c1, c2),
        }
    }
}

fn noncoplanar_circle_circle_intersect(c1: &Circle3d, c2: &Circle3d) -> CircleCircleIntersection {
    todo!()
}

/// Returns (larger, smaller)
/// If equal, returns (c1, c2)
fn cmp_radius_circle<'a>(c1: &'a Circle3d, c2: &'a Circle3d) -> (&'a Circle3d, &'a Circle3d) {
    match approx_cmp_f64(c1.radius(), c2.radius()) {
        FloatOrdering::Less => (c2, c1),
        _ => (c1, c2),
    }
}

fn coplanar_circle_circle_intersect(c1: &Circle3d, c2: &Circle3d) -> CircleCircleIntersection {
    let c1c2 = c2.center() - c1.center();
    let c1c2_norm_squared = c1c2.norm_squared();
    let r1r2_sum_squared = (c1.radius() + c2.radius()).powi(2);
    match approx_cmp_f64(c1c2_norm_squared, r1r2_sum_squared) {
        FloatOrdering::Equal => {
            let direction = UnitVector3::new_normalize(c1c2);
            let p = c1.center() + direction.scale(c1.radius());
            CircleCircleIntersection::Single(p)
        }
        FloatOrdering::Greater => CircleCircleIntersection::Empty,
        FloatOrdering::Less => {
            let r1r2_diff_squared = (c1.radius() - c2.radius()).powi(2);
            match approx_cmp_f64(c1c2_norm_squared, r1r2_diff_squared) {
                FloatOrdering::Less => CircleCircleIntersection::Empty,
                FloatOrdering::Equal => {
                    let (larger_c, smaller_c) = cmp_radius_circle(c1, c2);
                    let direction =
                        UnitVector3::new_normalize(smaller_c.center() - larger_c.center());
                    let p = larger_c.center() + direction.scale(larger_c.radius());
                    CircleCircleIntersection::Single(p)
                }
                FloatOrdering::Greater => {
                    let c1c2_normalized = UnitVector3::new_normalize(c1c2);
                    let c1c2_perpendicular = c1.n().cross(&c1c2_normalized);
                    // q = d^2 + r_1^2 - r_2^2
                    let q = c1c2_norm_squared + c1.radius().powi(2) - c2.radius().powi(2);
                    let c1c2_norm = c1c2.norm();
                    let dx = q / (2.0 * c1c2_norm);
                    let dy = (4.0 * c1c2_norm_squared * c1.radius().powi(2) - q).sqrt()
                        / (2.0 * c1c2_norm);
                    let p_dx = c1.center() + c1c2_normalized.scale(dx);
                    let p1 = p_dx + c1c2_perpendicular.scale(dy);
                    let p2 = p_dx - c1c2_perpendicular.scale(dy);
                    CircleCircleIntersection::Double(p1, p2)
                }
            }
        }
    }
}
