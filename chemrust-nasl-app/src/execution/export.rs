#[cfg(not(test))]
use std::fs::write;
use std::path::PathBuf;
use std::{fmt::Display, path::Path};

use castep_periodic_table::element::ElementSymbol;
use chemrust_core::data::lattice::CrystalModel;
use chemrust_nasl::{CoordSite, SearchReports, Visualize};

use crate::error::RunError;

#[derive(Debug, Clone)]
pub struct ExportFile<T: Display, P: AsRef<Path>> {
    file: T,
    file_name: P,
}

impl<T: Display, P: AsRef<Path>> ExportFile<T, P> {
    pub fn new(file: T, file_name: P) -> Self {
        Self { file, file_name }
    }

    pub fn file(&self) -> &T {
        &self.file
    }

    pub fn file_name(&self) -> &P {
        &self.file_name
    }
}

pub trait RhinoExport {
    fn new_element(&self) -> &ElementSymbol;
    fn base_model_name(&self) -> &str;
    fn export_dir(&self) -> &Path;
    fn potential_loc(&self) -> &Path;
    fn export_filename<E: ExportFormat>(&self, coord_site: &impl CoordSite) -> PathBuf {
        let atom_ids_text = coord_site.connecting_atoms_msg();
        self.export_dir().join(format!(
            "{}_{}.{}",
            self.base_model_name(),
            atom_ids_text,
            E::suffix()
        ))
    }
    fn create_all_sites<E: ExportFormat>(
        &self,
        base_model: &impl CrystalModel,
        coord_sites: &[(impl CoordSite + Visualize)],
    ) -> Vec<ExportFile<E::Item, PathBuf>> {
        coord_sites
            .iter()
            .map(|c| self.create_each_site::<E>(base_model, c))
            .collect()
    }
    fn create_each_site<E: ExportFormat>(
        &self,
        base_model: &impl CrystalModel,
        coord_site: &(impl CoordSite + Visualize),
    ) -> ExportFile<E::Item, PathBuf> {
        ExportFile::new(
            E::add_new_site(base_model, coord_site, self.new_element()),
            self.export_filename::<E>(coord_site),
        )
    }
    fn export_all_kinds_sites<E: ExportFormat>(
        &self,
        base_model: &impl CrystalModel,
        search_reports: &SearchReports,
    ) -> Result<(), RunError> {
        if let Some(points) = search_reports.viable_single_points() {
            self.create_all_sites::<E>(base_model, points)
                .iter()
                .try_for_each(|e| self.export_each_site::<E, PathBuf>(e))?;
        }
        if let Some(points) = search_reports.viable_double_points() {
            self.create_all_sites::<E>(base_model, points)
                .iter()
                .try_for_each(|e| self.export_each_site::<E, PathBuf>(e))?;
        }
        if let Some(points) = search_reports.points() {
            self.create_all_sites::<E>(base_model, points)
                .iter()
                .try_for_each(|e| self.export_each_site::<E, PathBuf>(e))?;
        }
        Ok(())
    }
    fn export_each_site<E: ExportFormat, P: AsRef<Path>>(
        &self,
        new_file: &ExportFile<E::Item, P>,
    ) -> Result<(), RunError> {
        #[cfg(test)]
        {
            println!("{}", new_file.file_name().as_ref().display());
            Ok(())
        }
        #[cfg(not(test))]
        {
            Ok(write(new_file.file_name(), new_file.file().to_string())?)
        }
    }
}

pub trait ExportFormat {
    type Item: Display;
    fn suffix() -> String;
    fn add_new_site(
        base_model: &impl CrystalModel,
        coord_site: &(impl CoordSite + Visualize),
        element_symbol: &ElementSymbol,
    ) -> Self::Item;
}

#[cfg(test)]
mod test {
    use crate::execution::search_job::SearchJob;
    use castep_cell_io::cell_document::CellDocument;
    use crystal_cif_io::DataBlock;

    use crate::{error::RunError, ModelFormat, TaskTable};

    use super::RhinoExport;

    #[test]
    fn export_traits() -> Result<(), RunError> {
        let table_path = "example_task.yaml";
        let task_table = TaskTable::load_task_table(table_path).expect("Path not found");
        let model = ModelFormat::load_model(task_table.model_path())
            .ok()
            .unwrap();
        let cell = model.as_cell().unwrap();
        let results = task_table.search_config().search(cell)?;
        assert!(task_table
            .export_all_kinds_sites::<CellDocument>(cell, &results)
            .is_ok());
        assert!(task_table
            .export_all_kinds_sites::<DataBlock>(cell, &results)
            .is_ok());
        Ok(())
    }
}
