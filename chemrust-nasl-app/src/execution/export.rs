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
        coord_sites: &[impl CoordSite + Visualize],
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

/// The CIF export item. The `ExportFormat` impl builds one
/// `crystal-cif-io` 0.2 data block and writes it as CIF text.
pub struct CifExport;

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
    use castep_cell_io::cell_document::CellDocument;
    use crystal_cif_io::parse_core_cif;

    use crate::execution::search_job::SearchJob;

    use crate::{
        error::{FormatError, RunError},
        execution::cif_model::CifModel,
        ModelFormat, TaskTable,
    };

    use super::{CifExport, ExportFormat, RhinoExport};

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
            .export_all_kinds_sites::<CifExport>(cell, &results)
            .is_ok());
        Ok(())
    }

    fn counts(results: &chemrust_nasl::SearchReports) -> (usize, usize, usize) {
        (
            results.viable_single_points().map(|v| v.len()).unwrap_or(0),
            results.viable_double_points().map(|v| v.len()).unwrap_or(0),
            results.points().map(|v| v.len()).unwrap_or(0),
        )
    }

    #[test]
    fn cif_model_matches_cell_model() -> Result<(), RunError> {
        let table = TaskTable::load_task_table("example_task.yaml")?;
        let config = table.search_config();
        let cif = ModelFormat::load_model("../scanner_test_models/CuTP001.cif")?
            .as_cif()
            .expect("the model loads as a CIF data block")
            .clone();
        let cell = ModelFormat::load_model("../scanner_test_models/CuTP001.cell")?
            .as_cell()
            .expect("the model loads as a CASTEP cell")
            .clone();
        let cif_results = config.search(&cif).ok();
        let cell_results = config.search(&cell).ok();
        let cif_counts = cif_results.as_ref().map(counts).unwrap_or((0, 0, 0));
        let cell_counts = cell_results.as_ref().map(counts).unwrap_or((0, 0, 0));
        dbg!(cif_counts, cell_counts);
        assert_eq!(cif_counts, cell_counts);
        Ok(())
    }

    #[test]
    fn cif_round_trip_search_matches_cell_search() -> Result<(), RunError> {
        let table = TaskTable::load_task_table("example_task.yaml")?;
        let config = table.search_config();
        let model = ModelFormat::load_model(table.model_path())?
            .as_cell()
            .expect("the example task points at a .cell model")
            .clone();
        let results = config.search(&model)?;
        let single_sites = results
            .viable_single_points()
            .ok_or(RunError::NoAvailableResults)?;
        let site = single_sites
            .first()
            .expect("the example task finds single sites");
        let cif_text = CifExport::add_new_site(&model, site, &config.new_element);
        let document = parse_core_cif(&cif_text)
            .map_err(|_| RunError::Message("the exported CIF does not parse".to_string()))?;
        let mut block = document
            .blocks
            .first()
            .cloned()
            .ok_or(RunError::FormatError(FormatError::Missing))?;
        // Drop the exported adsorbate row. The search sees the base model.
        block.atom_site.pop();
        let cif_model =
            CifModel::from_block(&block).map_err(|error| RunError::Message(error.to_string()))?;
        let cif_results = config.search(&cif_model)?;
        // SAC_GDY_V at 2.5 A bond length, z in [0.4, 0.5], finds
        // 24 single, 51 double, and 19 multi sites.
        assert_eq!(counts(&results), (24, 51, 19));
        assert_eq!(counts(&cif_results), counts(&results));
        Ok(())
    }
}
