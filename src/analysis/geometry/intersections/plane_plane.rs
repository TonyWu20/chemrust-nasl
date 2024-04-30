use std::f64::EPSILON;

use nalgebra::{Point3, UnitVector3};

use crate::analysis::geometry::primitives::{Line, Plane};

pub fn plane_plane_intersection(p1: &Plane, p2: &Plane) -> Option<Line> {
    let n1 = p1.normal();
    let n2 = p2.normal();
    let n3 = n1.cross(&n2);
    if n3.norm_squared().abs() < EPSILON {
        return None;
    }
    let h1 = -p1.d();
    let h2 = -p2.d();
    let d1 = (h1 - h2 * n1.dot(&n2)) / (1.0 - n1.dot(&n2).powi(2));
    let d2 = (h2 - h1 * n1.dot(&n2)) / (1.0 - n1.dot(&n2).powi(2));
    let x0 = Point3::from(n1.scale(d1) + n2.scale(d2));
    let line = Line::new(x0, UnitVector3::new_normalize(n3));
    Some(line)
}

#[cfg(test)]
mod circle_circle {
    use std::f64::EPSILON;

    use nalgebra::{Matrix3x2, Point3, UnitVector3, Vector2, Vector3};

    use crate::analysis::geometry::primitives::{Circle3d, Plane};

    use super::plane_plane_intersection;

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
        let p_p = plane_plane_intersection(&plane_1, &plane_2);
        println!("{:?}", p_p);
        let plane_3 = Plane::new(UnitVector3::new_normalize(Vector3::new(1.0, 0.0, 0.0)), 3.0);
        let plane_4 = Plane::new(UnitVector3::new_normalize(Vector3::new(0.0, 1.0, 0.0)), 6.0);
        let p_p = plane_plane_intersection(&plane_3, &plane_4);
        println!("{:?}", p_p);
        let plane_5 = Plane::new(UnitVector3::new_normalize(Vector3::new(1.0, 0.0, 0.0)), 3.0);
        let plane_6 = Plane::new(UnitVector3::new_normalize(Vector3::new(1.0, 0.0, 0.0)), 6.0);
        let p_p = plane_plane_intersection(&plane_5, &plane_6);
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
