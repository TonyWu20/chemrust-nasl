use std::{error::Error, path::Path};

use chemrust_core::data::{
    atom::CoreAtomData,
    geom::coordinates::CoordData,
    lattice::{cell_param::UnitCellParameters, CrystalModel, LatticeCell},
};
use nalgebra::Point3;

use crate::supportive_data::FractionalCoordRange;

use super::{
    format_identify::{self, match_format},
    format_loader::load_cell_file,
};

pub fn boundary_check(v: f64) -> f64 {
    if v < 0.0 {
        v + 1.0
    } else if v > 1.0 {
        v - 1.0
    } else {
        v
    }
}
pub fn load_model<P: AsRef<Path>>(model_path: &P) -> Result<LatticeCell, Box<dyn Error>> {
    let format = match_format(model_path)?;

    match format {
        format_identify::AcceptFormat::Cell => load_cell_file(model_path),
    }
}

pub fn get_to_check_atom(
    model: &LatticeCell,
    x_range: FractionalCoordRange,
    y_range: FractionalCoordRange,
    z_range: FractionalCoordRange,
) -> Vec<(usize, Point3<f64>)> {
    model
        .get_atom_data()
        .coords()
        .iter()
        .filter_map(|cd| match cd {
            CoordData::Fractional(frac) => {
                let point = frac.map(boundary_check);
                if x_range.is_in_range(point.x)
                    && y_range.is_in_range(point.y)
                    && z_range.is_in_range(point.z)
                {
                    Some(model.get_cell_parameters().cell_tensor() * point)
                } else {
                    None
                }
            }
            CoordData::Cartesian(cart) => {
                let frac = model
                    .get_cell_parameters()
                    .cell_tensor()
                    .try_inverse()
                    .unwrap()
                    * cart;
                let point = frac.map(boundary_check);
                if x_range.is_in_range(point.x)
                    && y_range.is_in_range(point.y)
                    && z_range.is_in_range(point.z)
                {
                    Some(model.get_cell_parameters().cell_tensor() * point)
                } else {
                    None
                }
            }
        })
        .enumerate()
        .collect()
}
