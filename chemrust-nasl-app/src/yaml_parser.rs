use std::{
    fs::read_to_string,
    path::{Path, PathBuf},
};

use castep_cell_io::{
    cell_document::CellDocument, CastepParams, CastepTask, CellParser, EnergyCutoffError,
};
use castep_periodic_table::{
    data::ELEMENT_TABLE,
    element::{Element, ElementSymbol, LookupElement},
};
use castep_seeding::{CellBuilding, ParamBuilding, RootJobs, SeedFolder, SeedingErrors};
use serde::{Deserialize, Serialize};

use crate::{
    error::RunError,
    execution::{RhinoExport, SearchJob},
    interactive_ui::KPointQuality,
    supportive_data::FractionalCoordRange,
};
pub use export_config::ExportConfig;
pub use run_config::RunConfig;

mod export_config;
mod run_config;

#[derive(Debug, Serialize, Deserialize, Clone)]
/// A config struct
pub struct TaskTable {
    #[serde(flatten)]
    pub(crate) search_config: RunConfig,
    #[serde(flatten)]
    pub(crate) export_config: ExportConfig,
}

impl TaskTable {
    pub fn new(search_config: RunConfig, export_config: ExportConfig) -> Self {
        Self {
            search_config,
            export_config,
        }
    }

    pub fn load_task_table<P: AsRef<Path>>(filepath: P) -> Result<Self, RunError> {
        let table_src = std::fs::File::open(filepath)?;
        let table = serde_yaml::from_reader(table_src)?;
        Ok(table)
    }

    pub fn model_path(&self) -> &str {
        self.search_config.model_path.as_ref()
    }

    pub fn new_element(&self) -> &Element {
        ELEMENT_TABLE.get_by_symbol(self.search_config.new_element)
    }

    pub fn target_bondlength(&self) -> f64 {
        self.search_config.target_bondlength
    }

    pub fn export_dir(&self) -> &PathBuf {
        &self.export_config.export_dir
    }

    pub fn kpoint_quality(&self) -> &KPointQuality {
        &self.export_config.kpoint_quality
    }

    pub fn edft(&self) -> bool {
        self.export_config.edft
    }

    pub fn potential_dir(&self) -> Option<&PathBuf> {
        self.export_config.potential_dir.as_ref()
    }
    pub fn x_range(&self) -> FractionalCoordRange {
        FractionalCoordRange::new(self.search_config.x_range.0, self.search_config.x_range.1)
    }
    pub fn y_range(&self) -> FractionalCoordRange {
        FractionalCoordRange::new(self.search_config.y_range.0, self.search_config.y_range.1)
    }
    pub fn z_range(&self) -> FractionalCoordRange {
        FractionalCoordRange::new(self.search_config.z_range.0, self.search_config.z_range.1)
    }

    pub fn search_config(&self) -> &RunConfig {
        &self.search_config
    }

    pub fn export_config(&self) -> &ExportConfig {
        &self.export_config
    }
}

impl SearchJob for RunConfig {
    fn target_bondlength(&self) -> f64 {
        self.target_bondlength
    }

    fn search_range(&self, axis: crate::execution::Axis) -> FractionalCoordRange {
        match axis {
            crate::execution::Axis::X => self.x_range(),
            crate::execution::Axis::Y => self.y_range(),
            crate::execution::Axis::Z => self.z_range(),
        }
    }
}

impl RhinoExport for TaskTable {
    fn new_element(&self) -> &ElementSymbol {
        &self.search_config.new_element
    }

    fn base_model_name(&self) -> &str {
        Path::new(self.model_path())
            .file_stem()
            .expect("Given path have a proper file stem")
            .to_str()
            .expect("Filename contains valid UTF-8 characters only")
    }

    fn export_dir(&self) -> &Path {
        self.export_dir()
    }

    /// Use local directories by default if not provided
    fn potential_loc(&self) -> &Path {
        self.potential_dir().map_or(Path::new("Potentials"), |v| v)
    }
}

struct Seed<P: AsRef<Path>> {
    cell_path: P,
    cell_doc: CellDocument,
}

impl<P: AsRef<Path>> Seed<P> {
    fn from_cell_path(cell_path: P) -> Result<Self, SeedingErrors> {
        let cell_content = read_to_string(&cell_path).map_err(SeedingErrors::ReadToString)?;
        let cell_doc = CellParser::from(&cell_content)
            .parse()
            .map_err(SeedingErrors::CellParseError)?;
        Ok(Self {
            cell_path,
            cell_doc,
        })
    }
}

impl<P> SeedFolder for Seed<P>
where
    P: AsRef<Path>,
{
    fn seed_name(&self) -> &str {
        self.cell_path
            .as_ref()
            .file_stem()
            .unwrap()
            .to_str()
            .unwrap()
    }

    fn root_dir(&self) -> impl AsRef<Path> {
        self.cell_path.as_ref().parent().unwrap()
    }

    fn cell_template(&self) -> &castep_cell_io::cell_document::CellDocument {
        &self.cell_doc
    }
}

impl RootJobs for TaskTable {
    fn root_path(&self) -> impl AsRef<Path> {
        self.export_dir()
    }

    fn generate_seed_folders(
        &self,
    ) -> Result<Vec<impl castep_seeding::SeedFolder>, castep_seeding::SeedingErrors> {
        self.get_cell_paths()?
            .into_iter()
            .map(Seed::from_cell_path)
            .collect()
    }
}

impl CellBuilding for ExportConfig {}

impl ParamBuilding for ExportConfig {
    fn build_param_for_task(
        &self,
        template_cell: &CellDocument,
        castep_task: CastepTask,
    ) -> Result<CastepParams, EnergyCutoffError> {
        match castep_task {
            castep_cell_io::CastepTask::BandStructure => self.bs_param_template(
                template_cell,
                castep_cell_io::EnergyCutoff::Ultrafine,
                self.edft,
                self.potential_loc(),
            ),
            castep_cell_io::CastepTask::GeometryOptimization => self.geom_opt_param_template(
                template_cell,
                castep_cell_io::EnergyCutoff::Ultrafine,
                self.edft,
                self.potential_loc(),
            ),
        }
    }
}

#[cfg(test)]
mod test {
    use super::TaskTable;

    #[test]
    fn test_task_table() {
        let table_path = "example_task.yaml";
        let task_table = TaskTable::load_task_table(table_path).expect("Path not found");
        println!("{}", task_table.model_path());
        println!("{}", task_table.kpoint_quality());
        println!("{:#?}", task_table.x_range());
        println!(
            "{}",
            task_table
                .export_dir()
                .file_name()
                .unwrap()
                .to_str()
                .unwrap()
        );
    }
}
