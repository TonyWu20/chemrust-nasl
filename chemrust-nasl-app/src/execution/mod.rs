use std::{error::Error, fs::read_to_string};

use castep_cell_io::CellParser;
use chemrust_core::data::{
    atom::CoreAtomData,
    geom::coordinates::CoordData,
    lattice::{cell_param::UnitCellParameters, CrystalModel},
};
use chemrust_nasl::{search_sites, SearchConfig, SearchReports, SiteIndex};
use nalgebra::Point3;

use crate::{error::RunError, yaml_parser::TaskTable};

use self::{
    export::export_all,
    format_loader::load_cell_file,
    helpers::{boundary_check, get_to_check_atom, load_model},
};

mod export;
mod format_identify;
mod format_loader;
mod helpers;

pub fn search(task_config: &TaskTable) -> Result<SearchReports, Box<dyn Error>> {
    let model = load_model(&task_config.model_path())?;
    let to_check = get_to_check_atom(
        &model,
        task_config.x_range(),
        task_config.y_range(),
        task_config.z_range(),
    );
    let points: Vec<Point3<f64>> = model
        .get_atom_data()
        .coords()
        .iter()
        .map(|cd| match cd {
            CoordData::Fractional(frac) => {
                let point = frac.map(boundary_check);
                model.get_cell_parameters().cell_tensor() * point
            }
            // Todo: boundary check for cartesian coordinates
            CoordData::Cartesian(p) => *p,
        })
        .collect();
    let site_index = SiteIndex::new(points);
    let search_config = SearchConfig::new(&to_check, task_config.target_bondlength());
    let search_report = search_sites(&site_index, &search_config);
    if search_report.viable_single_points().is_none()
        && search_report.viable_double_points().is_none()
        && search_report.points().is_none()
    {
        Err(Box::new(RunError::Message(
            "No available results for this config.".to_string(),
        )))
    } else {
        Ok(search_report)
    }
}

pub fn export_results_in_cell(
    task_config: &TaskTable,
    search_results: &SearchReports,
) -> Result<(), Box<dyn Error>> {
    let content = read_to_string(&task_config.model_path)?;
    let base_model = CellParser::from(content.as_str()).parse()?;
    let cell = load_cell_file(&task_config.model_path)?;
    let cell_param = cell.get_cell_parameters();
    export_all(&base_model, cell_param, task_config, search_results)?;
    Ok(())
}
