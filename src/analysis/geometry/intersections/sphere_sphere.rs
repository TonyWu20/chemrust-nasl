use std::f64::EPSILON;

use nalgebra::{Point3, UnitVector3};

use crate::analysis::geometry::primitives::{Circle3d, Sphere};

#[derive(Debug)]
/// Relationship between two spheres
enum SphereSphereRelationship<'a> {
    /// d > r1 + r2
    TooFarAway,
    /// d = 0, r1 = r2
    Overlaps,
    /// d + r_smaller < r_larger
    InsideOutOfReach,
    /// d + r_smaller = r_larger, returns larger and direction from smaller center to larger center
    InsideCut(&'a Sphere, UnitVector3<f64>),
    /// d = r1 + r2
    OutsideCut,
    /// r_large - r_smaller < d < r1 + r2
    Intersect,
}

impl<'a> SphereSphereRelationship<'a> {
    /// returns (larger, smaller)
    /// if the two spheres are considered to be identical, returns the first as the larger
    fn cmp_sphere(s1: &'a Sphere, s2: &'a Sphere) -> (&'a Sphere, &'a Sphere) {
        if s1.radius() - s2.radius() > EPSILON {
            (s1, s2)
        } else if s1.radius() - s2.radius() < -1.0 * EPSILON {
            (s2, s1)
        } else {
            (s1, s2)
        }
    }
    /// Determine the relationship
    fn determine(s1: &'a Sphere, s2: &'a Sphere) -> Self {
        let d = s2.center() - s1.center();
        let d_norm2 = d.norm_squared();
        let r1_plus_r2_2 = (s1.radius() + s2.radius()).powi(2);
        let r1_diff_r2_2 = (s1.radius() - s2.radius()).powi(2);
        let larger_r = f64::max(s1.radius(), s2.radius());
        let smaller_r = f64::min(s1.radius(), s2.radius());
        // edge cases
        // 1. two spheres too far away (r1 + r2) < d
        if d_norm2 - r1_plus_r2_2 > EPSILON {
            Self::TooFarAway
        }
        // 2. one sphere is inside another, and the inner one
        // doesn't touch the outer.
        // d + r_small < r_large
        else if d_norm2 - (larger_r - smaller_r).powi(2) < -1.0 * EPSILON {
            Self::InsideOutOfReach
        }
        // 3. Two spheres touch from outside
        // d = r1 + r2
        // Take the floating point inaccracy in account
        else if (d_norm2 - r1_plus_r2_2).abs() < EPSILON {
            Self::OutsideCut
        }
        // 4. Overlaps
        // d = 0, r1-r2 = 0
        else if d_norm2 < EPSILON && r1_diff_r2_2 < EPSILON {
            Self::Overlaps
        }
        // 5. One touches the outer from inside
        // d + r_small = r_large
        else if (d_norm2 - r1_diff_r2_2).abs() < EPSILON {
            let (larger, smaller) = Self::cmp_sphere(s1, s2);
            let direction = UnitVector3::new_normalize(smaller.center() - larger.center());
            Self::InsideCut(larger, direction)
        } else {
            Self::Intersect
        }
    }
}

#[derive(Debug)]
/// Result of Sphere-Sphere Intersection
pub enum SphereSphereResult {
    Empty,
    Point(Point3<f64>),
    Circle(Circle3d),
    Overlap(Sphere),
}

pub fn sphere_sphere_intersect(s1: &Sphere, s2: &Sphere) -> SphereSphereResult {
    // S1->S2
    let d = s2.center() - s1.center();
    let sphere_intersect_cases = SphereSphereRelationship::determine(s1, s2);
    match sphere_intersect_cases {
        SphereSphereRelationship::Intersect => {
            // d1 = d/2 + (r1^2 - r2^2)/2d
            let d1 = 0.5 * d.norm() + 0.5 * (s1.radius().powi(2) - s2.radius().powi(2)) / d.norm();
            // h^2 = r1^2 - d1^2 = (r1+d1)(r1-d1)
            let h = ((s1.radius() + d1) * (s1.radius() - d1)).sqrt();
            let norm = UnitVector3::new_normalize(d);
            let center = s1.center() + norm.scale(d1);
            let circle = Circle3d::new(center, h, norm);
            SphereSphereResult::Circle(circle)
        }
        SphereSphereRelationship::OutsideCut => {
            SphereSphereResult::Point(s1.point_at_surface(&UnitVector3::new_normalize(d)))
        }
        SphereSphereRelationship::TooFarAway => SphereSphereResult::Empty,
        SphereSphereRelationship::Overlaps => SphereSphereResult::Overlap(*s1),
        SphereSphereRelationship::InsideOutOfReach => SphereSphereResult::Empty,
        SphereSphereRelationship::InsideCut(larger, direction) => {
            SphereSphereResult::Point(larger.point_at_surface(&direction))
        }
    }
}

#[cfg(test)]
mod test {

    use std::f64::EPSILON;

    use nalgebra::Point3;

    use crate::analysis::geometry::{
        intersections::sphere_sphere::{sphere_sphere_intersect, SphereSphereRelationship},
        primitives::Sphere,
    };

    #[test]
    fn sphere_sphere() {
        let s1 = Sphere::new(Point3::origin(), 2.0); // at origin, r = 2
        let s2 = Sphere::new(Point3::new(1.0, 1.0, 0.0), 2.0); // normal intersects with s1
        let s3 = Sphere::new(Point3::new(1.0, 0.0, 0.0), 1.0); // inside touches s1 at (2.0,0.0,0.0)
        let s4 = Sphere::new(Point3::origin(), 2.0 + EPSILON); // overlaps with s1
        let s5 = Sphere::new(Point3::new(0.0, 4.0, 0.0), 2.0 + EPSILON); // outside touches s1 (2.0,0.0,0.0)
        let s6 = Sphere::new(Point3::new(4.0, 0.0, 0.0), 1.999 + EPSILON); // outside empty with s1
        let s7 = Sphere::new(Point3::new(0.0, 0.0, 0.0), 1.0); // inside empty with s1
        let s8 = Sphere::new(Point3::new(2.0, 2.0, 2.0), 3.0);
        let pairs = [s2, s3, s4, s5, s6, s7, s8];
        pairs.iter().enumerate().for_each(|(id, s)| {
            println!(
                "\nCase {} result: {:?}",
                id,
                sphere_sphere_intersect(&s1, s)
            );
            println!(
                "Relationship: {:?}\n",
                SphereSphereRelationship::determine(&s1, s)
            )
        });
    }
}
