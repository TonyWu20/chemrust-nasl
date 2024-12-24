use chemrust_core::data::{
    atom::CoreAtomData,
    geom::coordinates::CoordData,
    lattice::{CrystalModel, UnitCellParameters},
};
use helpers::{boundary_check, get_inbound_cartesian_coord};
use std::fs::read_to_string;

use castep_cell_io::CellParser;
use chemrust_nasl::{SearchConfig, SearchReports};
use nalgebra::Point3;

use crate::{
    error::{FormatError, RunError},
    supportive_data::FractionalCoordRange,
    yaml_parser::TaskTable,
};

pub use export::export_all;
pub use format_identify::match_format;
pub use format_loader::{load_cell_content, load_cell_file};
pub use helpers::get_to_check_atom;

use self::helpers::load_model;

mod export;
mod format_identify;
mod format_loader;
mod helpers;

#[derive(Debug, Clone, Copy, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub enum Axis {
    X,
    Y,
    Z,
}

pub trait SearchJob {
    fn target_bondlength(&self) -> f64;
    fn search_range(&self, axis: Axis) -> FractionalCoordRange;
    fn get_to_check_atom<C: CrystalModel>(&self, model: &C) -> Vec<(usize, Point3<f64>)> {
        let x_range = self.search_range(Axis::X);
        let y_range = self.search_range(Axis::Y);
        let z_range = self.search_range(Axis::Z);
        model
            .get_atom_data()
            .coords_repr()
            .iter()
            .enumerate()
            .filter_map(|(i, cd)| {
                get_inbound_cartesian_coord(
                    cd,
                    x_range,
                    y_range,
                    z_range,
                    model.get_cell_parameters(),
                )
                .map(|p| (i, p))
            })
            .collect()
    }
    fn all_atoms<C: CrystalModel>(&self, model: &C) -> Vec<Point3<f64>> {
        let all_range = FractionalCoordRange::new(0.0, 1.0);
        model
            .get_atom_data()
            .coords_repr()
            .iter()
            .filter_map(|cd| {
                get_inbound_cartesian_coord(
                    cd,
                    all_range,
                    all_range,
                    all_range,
                    model.get_cell_parameters(),
                )
            })
            .collect()
    }
    fn search<C: CrystalModel>(&self, model: &C) -> Result<SearchReports, RunError> {
        let to_check = self.get_to_check_atom(model);
        let all_atoms = self.all_atoms(model);
        let search_config = SearchConfig::new(&to_check, &all_atoms, self.bondlength());
        let search_report = search_config.search_sites();
        if search_report.viable_single_points().is_none()
            && search_report.viable_double_points().is_none()
            && search_report.points().is_none()
        {
            Err(RunError::Message(
                "No available results for this config.".to_string(),
            ))
        } else {
            Ok(search_report)
        }
    }
    fn export_results(&self)
}

pub fn search_with_length<T: CrystalModel>(
    model: &T,
    bondlength: f64,
    (x_low, x_high): (f64, f64),
    (y_low, y_high): (f64, f64),
    (z_low, z_high): (f64, f64),
) -> Result<SearchReports, RunError> {
    let x_range = (x_low, x_high);
    let y_range = (y_low, y_high);
    let z_range = (z_low, z_high);
    let to_check = get_to_check_atom(model, x_range, y_range, z_range);
    let all_range = (0.0, 1.0);
    let all_points: Vec<Point3<f64>> = get_to_check_atom(model, all_range, all_range, all_range)
        .iter()
        .map(|(_i, point)| *point)
        .collect();
    let search_config = SearchConfig::new(&to_check, &all_points, bondlength);
    let search_report = search_config.search_sites();
    if search_report.viable_single_points().is_none()
        && search_report.viable_double_points().is_none()
        && search_report.points().is_none()
    {
        Err(RunError::Message(
            "No available results for this config.".to_string(),
        ))
    } else {
        Ok(search_report)
    }
}

/// Run search with `TaskTable`
pub fn search_with_task_table(task_config: &TaskTable) -> Result<SearchReports, RunError> {
    let model = load_model(task_config.model_path())?;
    let bondlength = task_config.target_bondlength();
    let x_range = task_config.x_range;
    let y_range = task_config.y_range;
    let z_range = task_config.z_range;
    match model {
        format_identify::ModelFormat::Cell(cell) => {
            search_with_length(&cell, bondlength, x_range, y_range, z_range)
        }
        format_identify::ModelFormat::CifDataBlock(datablock) => {
            search_with_length(&datablock, bondlength, x_range, y_range, z_range)
        }
    }
}

pub fn export_results_in_cell(
    task_config: &TaskTable,
    search_results: &SearchReports,
) -> Result<(usize, usize, usize), RunError> {
    let content = read_to_string(&task_config.model_path)
        .map_err(|_| RunError::FormatError(FormatError::ReadToString))?;
    let base_model = CellParser::from(&content)
        .parse()
        .map_err(|_| RunError::FormatError(FormatError::Compatible))?;
    let cell = load_cell_file(&task_config.model_path).map_err(RunError::FormatError)?;
    let cell_param = cell.get_cell_parameters();
    let (mul, single, double) = export_all(&base_model, cell_param, task_config, search_results)
        .map_err(|_| RunError::IO)?;
    Ok((mul, single, double))
}
