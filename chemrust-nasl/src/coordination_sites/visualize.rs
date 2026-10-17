use std::f64::consts::FRAC_PI_2;

use elements_rs::Element;
use nalgebra::{Matrix3, Point3, Vector3};
use strum::VariantArray;

use crate::DelegatePoint;

use super::{CoordCircle, CoordSphere, MultiCoordPoint};

/// Map a coordination number to a visual element.
///
/// For `n` in `2..=103` this returns the element with atomic number
/// `n` (`Element::VARIANTS[n - 1]`). `n >= 104` falls back to
/// [`Element::Np`].
///
/// The `n = 0` and `n = 1` cases mirror the legacy
/// `castep-periodic-table` mapping: its `ElementSymbol` discriminants
/// were `H = 0` and `He = 2` (the gap at `1` made
/// `ElementSymbol::try_from(1u8)` fail and fall back to
/// `ElementSymbol::W`).
fn element_by_coordination_number(n: usize) -> Element {
    match n {
        0 => Element::H,
        1 => Element::W,
        n if n < 104 => Element::VARIANTS.get(n - 1).copied().unwrap_or(Element::W),
        _ => Element::Np,
    }
}

pub trait Visualize {
    type Output;
    fn determine_coord(&self) -> Point3<f64>;
    fn element_by_cn_number(&self) -> Element;
    fn draw_with_element(&self, element_symbol: Element) -> Self::Output;
    fn fractional_coord(&self, cell_tensor: Matrix3<f64>) -> Point3<f64> {
        cell_tensor.try_inverse().expect("Matrix is not invertible") * self.determine_coord()
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Atom {
    symbol: Element,
    coord: Point3<f64>,
}

impl Atom {
    /// `coord` is the Cartesian position in Angstrom.
    pub fn new(symbol: Element, coord: Point3<f64>) -> Self {
        Self { symbol, coord }
    }

    pub fn coord(&self) -> Point3<f64> {
        self.coord
    }

    pub fn symbol(&self) -> Element {
        self.symbol
    }
}

impl Visualize for CoordSphere {
    type Output = Atom;

    fn draw_with_element(&self, element_symbol: Element) -> Self::Output {
        let coord = self.determine_coord();
        Atom::new(element_symbol, coord)
    }

    fn determine_coord(&self) -> Point3<f64> {
        let center = self.sphere.center();
        let z_shift = Vector3::z_axis().scale(self.sphere.radius());
        center + z_shift
    }

    fn element_by_cn_number(&self) -> Element {
        Element::Xe
    }
}

impl Visualize for CoordCircle {
    type Output = Atom;

    fn draw_with_element(&self, element_symbol: Element) -> Self::Output {
        Atom::new(element_symbol, self.determine_coord())
    }

    fn determine_coord(&self) -> Point3<f64> {
        self.circle().get_point_on_circle(FRAC_PI_2)
    }

    fn element_by_cn_number(&self) -> Element {
        Element::Ne
    }
}

impl Visualize for MultiCoordPoint {
    type Output = Atom;

    fn draw_with_element(&self, element_symbol: Element) -> Self::Output {
        Atom::new(element_symbol, self.determine_coord())
    }

    fn determine_coord(&self) -> Point3<f64> {
        self.point
    }

    fn element_by_cn_number(&self) -> Element {
        element_by_coordination_number(self.atom_ids().len())
    }
}

impl<const N: usize> Visualize for DelegatePoint<N> {
    type Output = Atom;

    fn determine_coord(&self) -> Point3<f64> {
        self.point
    }

    fn element_by_cn_number(&self) -> Element {
        Element::Xe
    }

    fn draw_with_element(&self, element_symbol: Element) -> Self::Output {
        Atom::new(element_symbol, self.determine_coord())
    }
}
