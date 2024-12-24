use std::path::Path;

use chemrust_core::data::{
    atom::CoreAtomData,
    geom::coordinates::CoordData,
    lattice::{CrystalModel, UnitCellParameters},
};
use nalgebra::Point3;

use crate::{error::RunError, supportive_data::FractionalCoordRange};

use super::{format_identify::ModelFormat, format_loader::load_cif_file};
use super::{
    format_identify::{self, match_format},
    format_loader::load_cell_file,
};

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

pub fn load_model<P: AsRef<Path>>(model_path: P) -> Result<ModelFormat, RunError> {
    let format = match_format(&model_path).map_err(|e| RunError::FormatError(e))?;

    match format {
        format_identify::AcceptFormat::Cell => Ok(ModelFormat::Cell(
            load_cell_file(model_path).map_err(|e| RunError::FormatError(e))?,
        )),
        format_identify::AcceptFormat::Cif => Ok(ModelFormat::CifDataBlock(
            load_cif_file(model_path).map_err(|e| RunError::FormatError(e))?,
        )),
    }
}

pub fn get_to_check_atom<T: CrystalModel>(
    model: &T,
    (x_low, x_high): (f64, f64),
    (y_low, y_high): (f64, f64),
    (z_low, z_high): (f64, f64),
) -> Vec<(usize, Point3<f64>)> {
    let x_range = FractionalCoordRange::new(x_low, x_high);
    let y_range = FractionalCoordRange::new(y_low, y_high);
    let z_range = FractionalCoordRange::new(z_low, z_high);
    model
        .get_atom_data()
        .coords_repr()
        .iter()
        .enumerate()
        .filter_map(|(i, cd)| match cd {
            CoordData::Fractional(frac) => {
                let point = frac.map(boundary_check);
                if x_range.is_in_range(point.x)
                    && y_range.is_in_range(point.y)
                    && z_range.is_in_range(point.z)
                {
                    let point = model.get_cell_parameters().lattice_bases() * point;
                    Some((i, point))
                } else {
                    None
                }
            }
            CoordData::Cartesian(cart) => {
                let frac = model
                    .get_cell_parameters()
                    .lattice_bases()
                    .try_inverse()
                    .unwrap()
                    * cart;
                let point = frac.map(boundary_check);
                if x_range.is_in_range(point.x)
                    && y_range.is_in_range(point.y)
                    && z_range.is_in_range(point.z)
                {
                    let point = model.get_cell_parameters().lattice_bases() * point;
                    Some((i, point))
                } else {
                    None
                }
            }
        })
        .collect()
}
