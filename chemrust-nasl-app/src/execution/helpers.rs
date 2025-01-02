use chemrust_core::data::{geom::coordinates::CoordData, lattice::UnitCellParameters};
use nalgebra::Point3;

use crate::supportive_data::FractionalCoordRange;

pub fn boundary_check(v: f64) -> f64 {
    if !(0.0..=1.0).contains(&v) {
        // If v = -0.9, then returns -0.9 - (-1.0) = 0.1
        // If v = -3.8, then returns -3.8 - (-4.0) = 0.2
        // If v = 2.8 then returns 2.8 - 2.0 = 0.8
        v - v.floor()
    } else {
        v
    }
}

pub fn get_inbound_cartesian_coord(
    cd: &CoordData,
    x_range: FractionalCoordRange,
    y_range: FractionalCoordRange,
    z_range: FractionalCoordRange,
    lattice_param: &impl UnitCellParameters,
) -> Option<Point3<f64>> {
    match cd {
        CoordData::Fractional(frac) => {
            let point = frac.map(boundary_check);
            if point
                .iter()
                .zip([x_range, y_range, z_range])
                .all(|(&v, range)| range.is_in_range(v))
            {
                let point = lattice_param.lattice_bases() * point;
                Some(point)
            } else {
                None
            }
        }
        CoordData::Cartesian(cart) => {
            let frac = lattice_param.lattice_bases().try_inverse().unwrap() * cart;
            let point = frac.map(boundary_check);
            if point
                .iter()
                .zip([x_range, y_range, z_range])
                .all(|(&v, range)| range.is_in_range(v))
            {
                let point = lattice_param.lattice_bases() * point;
                Some(point)
            } else {
                None
            }
        }
    }
}
