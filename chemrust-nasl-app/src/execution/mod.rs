use chemrust_core::data::lattice::CrystalModel;
use std::fs::read_to_string;

use castep_cell_io::CellParser;
use chemrust_nasl::{search_sites, SearchConfig, SearchReports, SiteIndex};
use nalgebra::Point3;

use crate::{
    error::{FormatError, RunError},
    supportive_data::FractionalCoordRange,
    yaml_parser::TaskTable,
};

use self::{
    export::export_all, format_identify::match_format, format_loader::load_cell_file,
    helpers::get_to_check_atom,
};

mod export;
mod format_identify;
mod format_loader;
mod helpers;

pub fn search_with_length<T: CrystalModel>(
    model: &T,
    bondlength: f64,
    x_range: FractionalCoordRange,
    y_range: FractionalCoordRange,
    z_range: FractionalCoordRange,
) -> Result<SearchReports, RunError> {
    let to_check = get_to_check_atom(model, x_range, y_range, z_range);
    let all_range = FractionalCoordRange::new(0.0, 1.0);
    let all_points: Vec<Point3<f64>> = get_to_check_atom(model, all_range, all_range, all_range)
        .iter()
        .map(|(_i, point)| *point)
        .collect();
    let site_index = SiteIndex::new(&all_points);
    let search_config = SearchConfig::new(&to_check, bondlength);
    let search_report = search_sites(&site_index, &search_config);
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

pub fn search(task_config: &TaskTable) -> Result<SearchReports, RunError> {
    // let ModelFormat::Cell(model) = load_model(&task_config.model_path())?;
    let format = match_format(&task_config.model_path()).map_err(RunError::FormatError)?;
    let search_report = match format {
        format_identify::AcceptFormat::Cell => search_with_length(
            &load_cell_file(task_config.model_path()).map_err(RunError::FormatError)?,
            task_config.target_bondlength(),
            task_config.x_range(),
            task_config.y_range(),
            task_config.z_range(),
        )?,
    };
    Ok(search_report)
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
