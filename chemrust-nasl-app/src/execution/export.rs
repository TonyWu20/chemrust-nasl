use std::io::Error as IoError;
use std::path::Path;

use castep_cell_io::{CellDocument, IonicPosition};
use chemrust_core::data::lattice::cell_param::UnitCellParameters;
use chemrust_nasl::{CoordSite, SearchResults, Visualize};

use crate::yaml_parser::TaskTable;

pub fn export_all<T: UnitCellParameters>(
    base_model: &CellDocument,
    cell_param: &T,
    task_config: &TaskTable,
    results: &SearchResults,
) -> Result<(), IoError> {
    export(base_model, cell_param, task_config, results.spheres())?;
    export(base_model, cell_param, task_config, results.circles())?;
    export(base_model, cell_param, task_config, results.points())?;
    Ok(())
}

fn export<T: CoordSite + Visualize, U: UnitCellParameters>(
    base_model: &CellDocument,
    cell_param: &U,
    task_config: &TaskTable,
    coord_sites: &[T],
) -> Result<(), IoError> {
    coord_sites.iter().try_for_each(|site| {
        let atom_ids_text = site.connecting_atoms_msg();
        let model_name = Path::new(task_config.model_path())
            .file_stem()
            .and_then(|s| s.to_str())
            .expect("Invalid filename");
        let filename = Path::new(task_config.export_dir())
            .join(format!("{}_{}.cell", model_name, atom_ids_text));
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
