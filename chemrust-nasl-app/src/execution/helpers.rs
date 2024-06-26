use std::{error::Error, path::Path};

use chemrust_core::data::{
    atom::CoreAtomData,
    geom::coordinates::CoordData,
    lattice::{CrystalModel, UnitCellParameters},
};
use nalgebra::Point3;

use crate::supportive_data::FractionalCoordRange;

use super::{
    format_identify::{self, match_format, ModelFormat},
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
pub fn load_model<P: AsRef<Path>>(model_path: &P) -> Result<ModelFormat, Box<dyn Error>> {
    let format = match_format(model_path)?;

    match format {
        format_identify::AcceptFormat::Cell => Ok(ModelFormat::Cell(load_cell_file(model_path)?)),
    }
}

pub fn get_to_check_atom(
    model: &impl CrystalModel,
    x_range: FractionalCoordRange,
    y_range: FractionalCoordRange,
    z_range: FractionalCoordRange,
) -> Vec<(usize, Point3<f64>)> {
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
