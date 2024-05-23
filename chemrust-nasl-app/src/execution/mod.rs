use std::{error::Error, fs::read_to_string};

use castep_cell_io::CellParser;
use chemrust_core::data::{
    atom::CoreAtomData,
    geom::coordinates::CoordData,
    lattice::{cell_param::UnitCellParameters, CrystalModel},
};
use chemrust_nasl::{search_sites, SearchConfig, SearchReports, SearchResults, SiteIndex};
use nalgebra::Point3;

use crate::yaml_parser::TaskTable;

use self::{
    export::export_all,
    format_loader::load_cell_file,
    helpers::{boundary_check, get_to_check_atom, load_model},
};

mod export;
mod format_identify;
mod format_loader;
mod helpers;

pub fn search(task_config: &TaskTable) -> Result<SearchResults, Box<dyn Error>> {
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
    let search_sites = search_sites(&site_index, &search_config);
    #[cfg(debug_assertions)]
    {
        if let SearchResults::Found(report) = &search_sites {
            let validated_spheres =
                SearchReports::validated_results(report.spheres(), &site_index, &search_config);
            dbg!(validated_spheres.len(), report.spheres().len());
            let validated_circles =
                SearchReports::validated_results(report.circles(), &site_index, &search_config);
            dbg!(validated_circles.len(), report.circles().len());
            let validated_points =
                SearchReports::validated_results(report.points(), &site_index, &search_config);
            dbg!(validated_points.len(), report.points().len());
        }
    }
    Ok(search_sites)
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
