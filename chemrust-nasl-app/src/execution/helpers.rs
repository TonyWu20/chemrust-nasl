use chemrust_geometry::{FracCoord, LatticeVectors};
use nalgebra::Point3;

use crate::supportive_data::FractionalCoordRange;

/// Convert one fractional coordinate to the inbound Cartesian point.
///
/// The coordinate is wrapped into [0, 1). A point outside the
/// fractional ranges returns `None`.
pub fn get_inbound_cartesian_coord(
    frac: &FracCoord,
    x_range: FractionalCoordRange,
    y_range: FractionalCoordRange,
    z_range: FractionalCoordRange,
    lattice: &LatticeVectors,
) -> Option<Point3<f64>> {
    let mut wrapped = *frac;
    wrapped.wrap();
    if wrapped
        .iter()
        .zip([x_range, y_range, z_range])
        .all(|(v, range)| range.is_in_range(*v))
    {
        Some(lattice.tensor() * wrapped.0)
    } else {
        None
    }
}
