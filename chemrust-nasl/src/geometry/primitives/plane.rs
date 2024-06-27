use nalgebra::{Point3, UnitVector3, Vector3};

#[derive(Debug, Clone, Copy)]
/// n1*x + n2*y + n3*z + d = 0
/// # Note
/// The normal is a unit vector
/// So if two planes have different d and cross product of two norms is a 0 vector,
/// they must be parallel to each other.
pub struct Plane {
    normal: UnitVector3<f64>,
    d: f64,
}

impl Plane {
    pub fn new(normal: UnitVector3<f64>, d: f64) -> Self {
        Self { normal, d }
    }
    pub fn from_normal_and_point(normal: UnitVector3<f64>, point: Point3<f64>) -> Self {
        let op: Vector3<f64> = point - Point3::origin();
        let d = normal.dot(&op) * -1.0;
        Self { normal, d }
    }

    pub fn normal(&self) -> UnitVector3<f64> {
        self.normal
    }
    pub fn d(&self) -> f64 {
        self.d
    }
    pub fn point_in_plane(&self, point: Point3<f64>) -> bool {
        let op = point - Point3::origin();
        // Original f64::EPSILON seems to be too strict
        (self.normal.dot(&op) + self.d).abs() < f64::EPSILON * 5.0
    }
}

#[cfg(test)]
mod plane_test {
    use nalgebra::{Point3, UnitVector3, Vector3};

    use super::Plane;

    #[test]
    fn point_in_plane() {
        let normal = UnitVector3::new_normalize(Vector3::new(2.0, 3.0, 2.0));
        let p = Point3::new(1.0, 2.0, 3.0);
        let plane = Plane::from_normal_and_point(normal, p);
        let tp1 = Point3::new(2.0, 2.0, 3.0);
        let tp2 = Point3::new(3.0, 3.0, -0.5);
        let tp3 = Point3::new(9.0, -1.0, -0.5);
        let tp4 = Point3::new(5.0, 0.0, 2.0);
        assert!(!plane.point_in_plane(tp1));
        assert!(plane.point_in_plane(tp2));
        assert!(plane.point_in_plane(tp3));
        assert!(plane.point_in_plane(tp4));
    }
}
