use std::fmt::Debug;
use std::fs::{create_dir_all, write};
use std::ops::ControlFlow;
use std::path::{Path, PathBuf};
use std::{fs::create_dir, io::Error as IoError};

use castep_cell_io::{CellDocument, IonicPosition};
use chemrust_core::data::lattice::UnitCellParameters;
use chemrust_nasl::{CoordSite, DelegatePoint, MultiCoordPoint, SearchReports, Visualize};
use crystal_cif_io::to_cif_document;

use crate::yaml_parser::TaskTable;

pub fn export_all<T: UnitCellParameters>(
    base_model: &CellDocument,
    cell_param: &T,
    task_config: &TaskTable,
    results: &SearchReports,
) -> Result<(usize, usize, usize), IoError> {
    let mut num_multi = 0_usize;
    let mut num_single = 0_usize;
    let mut num_double = 0_usize;
    if let Some(multi_points) = results.points() {
        let boundary_checked: Vec<MultiCoordPoint> =
            points_boundary_check(multi_points, cell_param);
        if boundary_checked.len() > 1 {
            export(base_model, cell_param, task_config, &boundary_checked)?;
            collectively_export(base_model, cell_param, task_config, &boundary_checked)?;
            num_multi = boundary_checked.len();
            println!(
                "Exported {} multi-coordinated positions;",
                boundary_checked.len()
            );
        }
    }
    if let Some(single_points) = results.viable_single_points() {
        let boundary_checked: Vec<DelegatePoint<1>> =
            points_boundary_check(single_points, cell_param);
        if boundary_checked.len() > 1 {
            export(base_model, cell_param, task_config, &boundary_checked)?;
            collectively_export(base_model, cell_param, task_config, &boundary_checked)?;
            num_single = boundary_checked.len();
            println!(
                "Exported {} possible singly-coordinated positions;",
                boundary_checked.len()
            )
        }
    }
    if let Some(double_points) = results.viable_double_points() {
        let boundary_checked: Vec<DelegatePoint<2>> =
            points_boundary_check(double_points, cell_param);
        if boundary_checked.len() > 1 {
            export(base_model, cell_param, task_config, &boundary_checked)?;
            collectively_export(base_model, cell_param, task_config, &boundary_checked)?;
            num_double = boundary_checked.len();
            println!(
                "Exported {} possible doubly-coordinated positions;",
                boundary_checked.len()
            );
        }
    }
    Ok((num_multi, num_single, num_double))
}

fn points_boundary_check<T: Visualize + Clone, U: UnitCellParameters>(
    points: &[T],
    cell_param: &U,
) -> Vec<T> {
    points
        .iter()
        .filter(|cp| {
            let frac_coord = cp.fractional_coord(cell_param.lattice_bases());
            let check = frac_coord.iter().try_for_each(|&v| {
                if !(0.0..=1.0).contains(&v) {
                    ControlFlow::Break(())
                } else {
                    ControlFlow::Continue(())
                }
            });
            matches!(check, ControlFlow::Continue(()))
        })
        .cloned()
        .collect::<Vec<T>>()
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
        create_dir_all(export_dir_path)?;
    }
    coord_sites.iter().try_for_each(|site| {
        let filename = export_filename(site, task_config);
        let mut new_model = base_model.clone();
        let new_pos_coordinate = site.fractional_coord(cell_param.lattice_bases());
        let new_pos = IonicPosition::new(
            task_config.new_element().symbol(),
            new_pos_coordinate.into(),
            None,
        );
        new_model
            .model_description_mut()
            .ionic_pos_block_mut()
            .positions_mut()
            .push(new_pos);
        let cif_file = to_cif_document(&new_model, filename.file_stem().unwrap().to_str().unwrap());
        let cif_filename = filename.with_extension("cif");
        new_model.write_out(filename)?;
        write(cif_filename, cif_file.to_string())
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
        let new_pos_coordinate = site.fractional_coord(cell_param.lattice_bases());
        let symbol = site.element_by_cn_number();
        let new_pos = IonicPosition::new(symbol, new_pos_coordinate.into(), None);
        new_model
            .model_description_mut()
            .ionic_pos_block_mut()
            .positions_mut()
            .push(new_pos);
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
