
use nalgebra::{Point3, UnitVector3};

use crate::geometry::primitives::{Line, Plane};

use super::{approx_cmp_f64, FloatOrdering, Intersect};

#[derive(Debug, Clone, Copy)]
pub enum PlanePlaneIntersection {
    Same,
    Parallel,
    Intersect(Line),
}

impl Intersect for Plane {
    type Output = PlanePlaneIntersection;

    fn intersect(&self, rhs: &Self) -> Self::Output {
        let n1 = self.normal();
        let n2 = rhs.normal();
        let n3 = n1.cross(&n2);
        if n3.norm_squared() < f64::EPSILON {
            let d1 = self.d();
            let d2 = self.d();
            if let FloatOrdering::Equal = approx_cmp_f64(d1, d2) {
                PlanePlaneIntersection::Same
            } else {
                PlanePlaneIntersection::Parallel
            }
        } else {
            let h1 = -self.d();
            let h2 = -rhs.d();
            let d1 = (h1 - h2 * n1.dot(&n2)) / (1.0 - n1.dot(&n2).powi(2));
            let d2 = (h2 - h1 * n1.dot(&n2)) / (1.0 - n1.dot(&n2).powi(2));
            let x0 = Point3::from(n1.scale(d1) + n2.scale(d2));
            let line = Line::new(x0, UnitVector3::new_normalize(n3));
            PlanePlaneIntersection::Intersect(line)
        }
    }
}

#[cfg(test)]
mod test {

    use nalgebra::{Point3, UnitVector3, Vector3};

    use crate::geometry::{
        intersections::{plane_plane::PlanePlaneIntersection, Intersect},
        primitives::{Circle3d, Plane},
    };

    #[test]
    fn plane_plane() {
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
        let plane_1 = c1.plane_of_circle();
        let plane_2 = c2.plane_of_circle();
        if let PlanePlaneIntersection::Intersect(line) = plane_1.intersect(&plane_2) {
            let origin = line.origin();
            let d = origin - c1.center();
            let ld = line.direction();
            let angle = d.angle(&ld);
            println!("{:?}", line);
            println!("{}", d);
            println!("{}", angle);
            let h = d.norm() * angle.sin();
            println!("{}, {}", h, d.norm());
        }
        let plane_3 = Plane::new(UnitVector3::new_normalize(Vector3::new(1.0, 0.0, 0.0)), 3.0);
        let plane_4 = Plane::new(UnitVector3::new_normalize(Vector3::new(0.0, 1.0, 0.0)), 6.0);
        let p_p = plane_3.intersect(&plane_4);
        println!("{:?}", p_p);
        let plane_5 = Plane::new(UnitVector3::new_normalize(Vector3::new(1.0, 0.0, 0.0)), 3.0);
        let plane_6 = Plane::new(UnitVector3::new_normalize(Vector3::new(1.0, 0.0, 0.0)), 6.0);
        let p_p = plane_5.intersect(&plane_6);
        println!("{:?}", p_p);
    }
    #[test]
    fn coplanar() {
        let n1 = Vector3::new(1.0, 0.0, 0.0);
        let n2 = Vector3::new(1.5, 0.0, 0.0);
        let n3 = n1.cross(&n2);
        println!("{}", n3.norm_squared());
    }
}
