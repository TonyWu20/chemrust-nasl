use nalgebra::Point3;

use crate::geometry::{
    approx_cmp_f64, intersections::circle_circle::coplanar_circle_circle_intersect, Circle3d,
    CircleCircleIntersection, Sphere,
};

use super::{FloatOrdering, Intersect};

#[derive(Debug, Clone, Copy)]
pub enum CircleSphereIntersection {
    Invalid,
    Zero,
    InsideSphere,
    SphereInCircle,
    Single(Point3<f64>),
    Double(Point3<f64>, Point3<f64>),
    Circle(Circle3d),
}

impl Intersect<Sphere> for Circle3d {
    type Output = CircleSphereIntersection;

    fn intersect(&self, rhs: &Sphere) -> Self::Output {
        let cs_cc = self.center() - rhs.center();
        // the n is unit vector so this means the projection distance of vector
        // cs_cc on the unit normal vector of the circle plane.
        let cut_at = self.n().dot(&cs_cc);
        // If the absolute value of projection distance is greater than the sphere
        // radius, the intersection plane is above or below the sphere. No
        // intersection.
        match approx_cmp_f64(cut_at.abs(), rhs.radius()) {
            FloatOrdering::Greater => CircleSphereIntersection::Zero,
            FloatOrdering::Equal => {
                let projected_dist = (cs_cc.norm_squared() - rhs.radius().powi(2)).sqrt();
                match approx_cmp_f64(projected_dist, self.radius()) {
                    FloatOrdering::Less => CircleSphereIntersection::SphereInCircle,
                    FloatOrdering::Equal => {
                        let p = rhs.center() + self.n().scale(cut_at);
                        CircleSphereIntersection::Single(p)
                    }
                    FloatOrdering::Greater => CircleSphereIntersection::Zero,
                }
            }
            FloatOrdering::Less => {
                let new_circle_center = rhs.center() + self.n().scale(cut_at);
                // new circle radius <= Sphere radius
                let new_circle_radius = (rhs.radius().powi(2) - cut_at.powi(2)).sqrt();
                let new_circle = Circle3d::new(new_circle_center, new_circle_radius, self.n());
                let result = coplanar_circle_circle_intersect(self, &new_circle);
                match result {
                    CircleCircleIntersection::Empty => CircleSphereIntersection::Zero,
                    CircleCircleIntersection::Single(p) => CircleSphereIntersection::Single(p),
                    CircleCircleIntersection::Double(p1, p2) => {
                        CircleSphereIntersection::Double(p1, p2)
                    }
                    CircleCircleIntersection::Overlap(c) => CircleSphereIntersection::Circle(c),
                    CircleCircleIntersection::Contains(inner, _outer) => {
                        // No more float point calculation so exact number will
                        // be preserved
                        if inner.radius() == self.radius() {
                            CircleSphereIntersection::InsideSphere
                        } else {
                            CircleSphereIntersection::SphereInCircle
                        }
                    }
                    CircleCircleIntersection::Invalid => CircleSphereIntersection::Invalid,
                }
            }
        }
    }
}
