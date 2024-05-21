use castep_periodic_table::element::ElementSymbol;
use chemrust_core::data::geom::coordinates::CoordData;
use nalgebra::{UnitVector3, Vector3};

use super::{CoordCircle, CoordPoint, CoordSphere};

pub trait Visualize {
    type Output;
    fn draw_with_element(&self, element_symbol: &ElementSymbol) -> Self::Output;
}

pub struct Atom {
    symbol: ElementSymbol,
    coord: CoordData,
}

impl Atom {
    pub fn new(symbol: ElementSymbol, coord: CoordData) -> Self {
        Self { symbol, coord }
    }

    pub fn coord(&self) -> CoordData {
        self.coord
    }

    pub fn symbol(&self) -> ElementSymbol {
        self.symbol
    }
}

impl Visualize for CoordSphere {
    type Output = Atom;

    fn draw_with_element(&self, element_symbol: &ElementSymbol) -> Self::Output {
        let center = self.sphere.center();
        let z_shift = Vector3::z_axis().scale(self.sphere.radius());
        let coord = CoordData::Cartesian(center + z_shift);
        Atom::new(*element_symbol, coord)
    }
}

impl Visualize for CoordCircle {
    type Output = Atom;

    fn draw_with_element(&self, element_symbol: &ElementSymbol) -> Self::Output {
        let (x, y, _z) = (self.circle.n().x, self.circle.n().y, self.circle.n().z);
        // We want the v2 to act as the "z-axis" after transformation
        let pre_v1 = Vector3::new(-1.0 * y, x, 0.0);
        let v1 = if pre_v1.norm_squared() < f64::EPSILON {
            // such v1 becomes a null vector when the normal is (0, 0, z)
            Vector3::x_axis()
        } else {
            UnitVector3::new_normalize(pre_v1)
        };
        let v2 = self.circle.n().cross(&v1);
        let p = self.circle.center() + (v1.scale(0.0) + v2.scale(1.0)).scale(self.circle.radius());
        Atom::new(*element_symbol, CoordData::Cartesian(p))
    }
}

impl Visualize for CoordPoint {
    type Output = Atom;

    fn draw_with_element(&self, element_symbol: &ElementSymbol) -> Self::Output {
        Atom::new(*element_symbol, CoordData::Cartesian(self.point))
    }
}
