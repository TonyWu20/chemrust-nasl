use std::f64::consts::FRAC_PI_2;

use castep_periodic_table::element::ElementSymbol;
use chemrust_core::data::geom::coordinates::CoordData;
use nalgebra::{Matrix3, Point3, UnitVector3, Vector3};

use crate::DelegatePoint;

use super::{CoordCircle, CoordSphere, MultiCoordPoint};

pub trait Visualize {
    type Output;
    fn determine_coord(&self) -> Point3<f64>;
    fn element_by_cn_number(&self) -> ElementSymbol;
    fn draw_with_element(&self, element_symbol: ElementSymbol) -> Self::Output;
    fn fractional_coord(&self, cell_tensor: Matrix3<f64>) -> Point3<f64> {
        cell_tensor.try_inverse().expect("Matrix is not invertible") * self.determine_coord()
    }
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

    fn draw_with_element(&self, element_symbol: ElementSymbol) -> Self::Output {
        let coord = CoordData::Cartesian(self.determine_coord());
        Atom::new(element_symbol, coord)
    }

    fn determine_coord(&self) -> Point3<f64> {
        let center = self.sphere.center();
        let z_shift = Vector3::z_axis().scale(self.sphere.radius());
        center + z_shift
    }

    fn element_by_cn_number(&self) -> ElementSymbol {
        ElementSymbol::Xe
    }
}

impl Visualize for CoordCircle {
    type Output = Atom;

    fn draw_with_element(&self, element_symbol: ElementSymbol) -> Self::Output {
        Atom::new(element_symbol, CoordData::Cartesian(self.determine_coord()))
    }

    fn determine_coord(&self) -> Point3<f64> {
        self.circle().get_point_on_circle(FRAC_PI_2)
    }

    fn element_by_cn_number(&self) -> ElementSymbol {
        ElementSymbol::Ne
    }
}

impl Visualize for MultiCoordPoint {
    type Output = Atom;

    fn draw_with_element(&self, element_symbol: ElementSymbol) -> Self::Output {
        Atom::new(element_symbol, CoordData::Cartesian(self.determine_coord()))
    }

    fn determine_coord(&self) -> Point3<f64> {
        self.point
    }

    fn element_by_cn_number(&self) -> ElementSymbol {
        if self.atom_ids().len() < 104 {
            ElementSymbol::try_from(self.atom_ids().len() as u8).unwrap_or(ElementSymbol::W)
        } else {
            ElementSymbol::Np
        }
    }
}

impl<const N: usize> Visualize for DelegatePoint<N> {
    type Output = Atom;

    fn determine_coord(&self) -> Point3<f64> {
        self.point
    }

    fn element_by_cn_number(&self) -> ElementSymbol {
        ElementSymbol::Xe
    }

    fn draw_with_element(&self, element_symbol: ElementSymbol) -> Self::Output {
        Atom::new(element_symbol, CoordData::Cartesian(self.determine_coord()))
    }
}
