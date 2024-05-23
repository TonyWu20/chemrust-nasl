use std::fmt::Debug;
use std::ops::ControlFlow;
use std::path::{Path, PathBuf};
use std::{fs::create_dir, io::Error as IoError};

use castep_cell_io::{CellDocument, IonicPosition};
use chemrust_core::data::lattice::cell_param::UnitCellParameters;
use chemrust_nasl::{CoordPoint, CoordSite, SearchReports, Visualize};

use crate::yaml_parser::TaskTable;

pub fn export_all<T: UnitCellParameters>(
    base_model: &CellDocument,
    cell_param: &T,
    task_config: &TaskTable,
    results: &SearchReports,
) -> Result<(), IoError> {
    if !results.spheres().is_empty() {
        export(base_model, cell_param, task_config, results.spheres())?;
        collectively_export(base_model, cell_param, task_config, results.spheres())?;
    }
    if !results.circles().is_empty() {
        export(base_model, cell_param, task_config, results.circles())?;
        collectively_export(base_model, cell_param, task_config, results.circles())?;
    }
    if !results.points().is_empty() {
        let boundary_checked: Vec<CoordPoint> = results
            .points()
            .iter()
            .filter(|cp| {
                let frac_coord = cp.fractional_coord(cell_param.cell_tensor());
                let check = frac_coord.iter().try_for_each(|&v| {
                    if !(0.0..=1.0).contains(&v) {
                        ControlFlow::Break(v)
                    } else {
                        ControlFlow::Continue(())
                    }
                });
                matches!(check, ControlFlow::Continue(()))
            })
            .cloned()
            .collect();
        export(base_model, cell_param, task_config, &boundary_checked)?;
        collectively_export(base_model, cell_param, task_config, &boundary_checked)?;
    }
    Ok(())
}

fn export_filename<T: CoordSite>(coord_site: &T, task_config: &TaskTable) -> PathBuf {
    let atom_ids_text = coord_site.connecting_atoms_msg();
    let model_name = Path::new(task_config.model_path())
        .file_stem()
        .and_then(|s| s.to_str())
        .expect("Invalid filename");
    Path::new(task_config.export_dir()).join(format!("{}_{}.cell", model_name, atom_ids_text))
}

fn export<T: CoordSite + Visualize, U: UnitCellParameters>(
    base_model: &CellDocument,
    cell_param: &U,
    task_config: &TaskTable,
    coord_sites: &[T],
) -> Result<(), IoError> {
    let export_dir_path = Path::new(task_config.export_dir());
    if !export_dir_path.exists() {
        create_dir(export_dir_path)?;
    }
    coord_sites.iter().try_for_each(|site| {
        let filename = export_filename(site, task_config);
        let mut new_model = base_model.clone();
        let new_pos_coordinate = site.fractional_coord(cell_param.cell_tensor());
        let new_pos = IonicPosition::new(
            task_config.new_element().symbol(),
            new_pos_coordinate.into(),
            None,
        );
        new_model.append_position(new_pos);
        new_model.write_out(filename)
    })
}

fn collectively_export<T: CoordSite + Visualize + Debug, U: UnitCellParameters>(
    base_model: &CellDocument,
    cell_param: &U,
    task_config: &TaskTable,
    coord_sites: &[T],
) -> Result<(), IoError> {
    let export_dir_path = Path::new(task_config.export_dir());
    if !export_dir_path.exists() {
        create_dir(export_dir_path)?;
    }
    let mut new_model = base_model.clone();
    coord_sites.iter().for_each(|site| {
        let new_pos_coordinate = site.fractional_coord(cell_param.cell_tensor());
        let symbol = site.element_by_cn_number();
        let new_pos = IonicPosition::new(symbol, new_pos_coordinate.into(), None);
        new_model.append_position(new_pos);
    });
    let model_name = Path::new(task_config.model_path())
        .file_stem()
        .and_then(|s| s.to_str())
        .expect("Invalid filename");
    let filename = Path::new(task_config.export_dir()).join(format!(
        "{}_{}_all.cell",
        model_name,
        coord_sites[0].site_type()
    ));
    new_model.write_out(filename)
}
