use nalgebra::{Point3, UnitVector3};

use crate::geometry::{
    intersections::{approx_cmp_f64, FloatOrdering},
    primitives::{Circle3d, Line},
};

use super::{approx_eq_point_f64, plane_plane::PlanePlaneIntersection, FloatEq, Intersect};

#[derive(Debug, Clone, Copy)]
enum CircleCircleRelationship {
    Coplanar,
    ParallelPlane,
    PlaneIntersect(Line),
}

impl CircleCircleRelationship {
    fn determine(c1: &Circle3d, c2: &Circle3d) -> Self {
        let p1 = c1.plane_of_circle();
        let p2 = c2.plane_of_circle();
        match p1.intersect(&p2) {
            PlanePlaneIntersection::Same => Self::Coplanar,
            PlanePlaneIntersection::Parallel => Self::ParallelPlane,
            PlanePlaneIntersection::Intersect(line) => Self::PlaneIntersect(line),
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub enum CircleCircleIntersection {
    Empty,
    Single(Point3<f64>),
    Double(Point3<f64>, Point3<f64>),
    Overlap(Circle3d),
    /// inner and outer
    Contains(Circle3d, Circle3d),
}

#[derive(Debug, Clone, Copy)]
pub struct CircleCoplanarLine(Line);

#[derive(Debug, Clone, Copy)]
pub enum CircleLineIntersection {
    Empty,
    Single(Point3<f64>),
    Double(Point3<f64>, Point3<f64>),
}

impl Intersect for Circle3d {
    type Output = CircleCircleIntersection;

    fn intersect(&self, rhs: &Self) -> Self::Output {
        let case = CircleCircleRelationship::determine(self, rhs);
        match case {
            CircleCircleRelationship::ParallelPlane => CircleCircleIntersection::Empty,
            CircleCircleRelationship::Coplanar => coplanar_circle_circle_intersect(self, rhs),
            CircleCircleRelationship::PlaneIntersect(line) => {
                noncoplanar_circle_circle_intersect(self, rhs, &line)
            }
        }
    }
}

impl Intersect<CircleCoplanarLine> for Circle3d {
    type Output = CircleLineIntersection;

    fn intersect(&self, rhs: &CircleCoplanarLine) -> Self::Output {
        let line = rhs.0;
        let distance_to_line = line.point_to_line_distance(&self.center());
        match approx_cmp_f64(distance_to_line, self.radius()) {
            FloatOrdering::Less => {
                let line_origin_to_p = self.center() - line.origin();
                let angle = line.direction().angle(&line_origin_to_p);
                let h = (self.radius().powi(2) - distance_to_line.powi(2)).sqrt();
                let t1 = line_origin_to_p.norm() * angle.cos() + h;
                let t2 = line_origin_to_p.norm() * angle.cos() - h;
                let p1 = line.point(t1);
                let p2 = line.point(t2);
                CircleLineIntersection::Double(p1, p2)
            }
            FloatOrdering::Equal => CircleLineIntersection::Single(line.origin()),
            FloatOrdering::Greater => CircleLineIntersection::Empty,
        }
    }
}

fn noncoplanar_circle_circle_intersect(
    c1: &Circle3d,
    c2: &Circle3d,
    intersect_line: &Line,
) -> CircleCircleIntersection {
    let coplanar_line = CircleCoplanarLine(*intersect_line);
    let c1_line_result: CircleLineIntersection = c1.intersect(&coplanar_line);
    let c2_line_result: CircleLineIntersection = c2.intersect(&coplanar_line);
    match (c1_line_result, c2_line_result) {
        (CircleLineIntersection::Empty, _) => CircleCircleIntersection::Empty,
        (_, CircleLineIntersection::Empty) => CircleCircleIntersection::Empty,
        (CircleLineIntersection::Single(p1), CircleLineIntersection::Single(p2)) => {
            match approx_eq_point_f64(p1, p2) {
                FloatEq::Eq => CircleCircleIntersection::Single(p1),
                FloatEq::NotEq => CircleCircleIntersection::Empty,
            }
        }
        (CircleLineIntersection::Single(p1), CircleLineIntersection::Double(p2, p3)) => {
            if let FloatEq::Eq = approx_eq_point_f64(p1, p2) {
                CircleCircleIntersection::Single(p1)
            } else if let FloatEq::Eq = approx_eq_point_f64(p1, p3) {
                CircleCircleIntersection::Single(p1)
            } else {
                CircleCircleIntersection::Empty
            }
        }
        (CircleLineIntersection::Double(p2, p3), CircleLineIntersection::Single(p1)) => {
            if let FloatEq::Eq = approx_eq_point_f64(p1, p2) {
                CircleCircleIntersection::Single(p1)
            } else if let FloatEq::Eq = approx_eq_point_f64(p1, p3) {
                CircleCircleIntersection::Single(p1)
            } else {
                CircleCircleIntersection::Empty
            }
        }
        (CircleLineIntersection::Double(p1, p2), CircleLineIntersection::Double(p3, p4)) => {
            // Full cases discussion
            let case_1 = (approx_eq_point_f64(p1, p3), approx_eq_point_f64(p2, p4));
            // If the corresponding relationship is the other way.
            match case_1 {
                // p1 != p3, p2 == p4, then p2
                (FloatEq::NotEq, FloatEq::Eq) => CircleCircleIntersection::Single(p2),
                // p1 == p3, p2 != p4, then p1
                (FloatEq::Eq, FloatEq::NotEq) => CircleCircleIntersection::Single(p1),
                // p1 == p3 && p2 == p4, double
                (FloatEq::Eq, FloatEq::Eq) => CircleCircleIntersection::Double(p1, p2),
                // p1 != p3 && p2 != p4, could be mismatch, further discuss as case_2
                (FloatEq::NotEq, FloatEq::NotEq) => {
                    // Not inline for tidyness in reading
                    let case_2 = (approx_eq_point_f64(p1, p4), approx_eq_point_f64(p2, p3));
                    match case_2 {
                        // Like above 3 cases
                        (FloatEq::NotEq, FloatEq::Eq) => CircleCircleIntersection::Single(p2),
                        (FloatEq::Eq, FloatEq::NotEq) => CircleCircleIntersection::Single(p1),
                        (FloatEq::Eq, FloatEq::Eq) => CircleCircleIntersection::Double(p1, p2),
                        // (p1 != p3 && p2 != p4) && (p1 != p4 && p2 != p3),
                        // they don't have common points for sure now
                        (FloatEq::NotEq, FloatEq::NotEq) => CircleCircleIntersection::Empty,
                    }
                }
            }
        }
    }
}

/// Returns (larger, smaller)
/// If equal, returns (c1, c2)
fn cmp_radius_circle<'a>(c1: &'a Circle3d, c2: &'a Circle3d) -> (&'a Circle3d, &'a Circle3d) {
    match approx_cmp_f64(c1.radius(), c2.radius()) {
        FloatOrdering::Less => (c2, c1),
        _ => (c1, c2),
    }
}

/// Trivial math deduction
pub(crate) fn coplanar_circle_circle_intersect(
    c1: &Circle3d,
    c2: &Circle3d,
) -> CircleCircleIntersection {
    let c1c2 = c2.center() - c1.center();
    let c1c2_norm_squared = c1c2.norm_squared();
    // When the two circles have identical centers
    if let FloatOrdering::Equal = approx_cmp_f64(c1c2_norm_squared, 0.0) {
        match approx_cmp_f64(c1.radius(), c2.radius()) {
            FloatOrdering::Equal => CircleCircleIntersection::Overlap(*c1),
            FloatOrdering::Less => CircleCircleIntersection::Contains(*c1, *c2),
            FloatOrdering::Greater => CircleCircleIntersection::Contains(*c2, *c1),
        }
    }
    // The two circles' centers are separated
    else {
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
                        let c1c2_perpendicular =
                            UnitVector3::new_normalize(c1.n().cross(&c1c2_normalized));
                        // q = d^2 + r_1^2 - r_2^2
                        // let c1c2_norm = c1c2.norm();
                        let h = (c1c2_norm_squared + c1.radius().powi(2) - c2.radius().powi(2))
                            / (2.0 * c1c2.norm());
                        // let dy = (4.0 * c1c2_norm_squared * c1.radius().powi(2) - q).sqrt()
                        //     / (2.0 * c1c2_norm);
                        let dy = (c1.radius().powi(2) - h.powi(2)).sqrt();
                        let p_dx = c1.center() + c1c2_normalized.scale(h);
                        let p1 = p_dx + c1c2_perpendicular.scale(dy);
                        let p2 = p_dx - c1c2_perpendicular.scale(dy);
                        CircleCircleIntersection::Double(p1, p2)
                    }
                }
            }
        }
    }
}

#[cfg(test)]
/// Todo: write tests to cover all possible cases of circle-circle intersection
mod test {
    use nalgebra::{Point3, UnitVector3, Vector3};

    use crate::geometry::{intersections::Intersect, primitives::Circle3d};

    #[test]
    fn circle_circle() {
        let c1 = Circle3d::new(
            Point3::new(0.5833333333333333, 0.5833333333333333, 0.5833333333333333),
            1.726026264767,
            UnitVector3::new_normalize(Vector3::new(1.0, 1.0, 1.0)),
        );
        let c2 = Circle3d::new(
            Point3::new(0.5, 0.5, 0.0),
            1.8708286933869707,
            UnitVector3::new_normalize(Vector3::new(1.0, 1.0, 0.0)),
        );
        let res = c1.intersect(&c2);
        println!("{:?}", res);
    }
}
