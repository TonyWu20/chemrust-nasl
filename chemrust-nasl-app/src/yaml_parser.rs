use std::{
    error::Error,
    path::{Path, PathBuf},
};

use castep_periodic_table::{
    data::ELEMENT_TABLE,
    element::{Element, ElementSymbol, LookupElement},
};
use serde::{Deserialize, Serialize};

use crate::{interactive_ui::KPointQuality, supportive_data::FractionalCoordRange};

#[derive(Debug, Serialize, Deserialize, Clone)]
/// A config struct
pub struct TaskTable {
    pub(crate) model_path: String,
    pub(crate) new_element: ElementSymbol,
    pub(crate) target_bondlength: f64,
    pub(crate) x_range: (f64, f64),
    pub(crate) y_range: (f64, f64),
    pub(crate) z_range: (f64, f64),
    pub(crate) export_dir: PathBuf,
    pub(crate) potential_dir: Option<String>,
    pub(crate) kpoint_quality: KPointQuality,
    pub(crate) edft: bool,
}

impl TaskTable {
    pub fn new(
        model_path: String,
        new_element: ElementSymbol,
        target_bondlength: f64,
        x_range: (f64, f64),
        y_range: (f64, f64),
        z_range: (f64, f64),
        export_dir: PathBuf,
        potential_dir: Option<String>,
        kpoint_quality: KPointQuality,
        edft: bool,
    ) -> Self {
        Self {
            model_path,
            new_element,
            target_bondlength,
            x_range,
            y_range,
            z_range,
            export_dir,
            potential_dir,
            kpoint_quality,
            edft,
        }
    }

    pub fn load_task_table<P: AsRef<Path>>(filepath: P) -> Result<Self, Box<dyn Error>> {
        let table_src = std::fs::File::open(filepath)?;
        let table = serde_yaml::from_reader(table_src)?;
        Ok(table)
    }

    pub fn model_path(&self) -> &str {
        self.model_path.as_ref()
    }

    pub fn new_element(&self) -> &Element {
        ELEMENT_TABLE.get_by_symbol(self.new_element)
    }

    pub fn target_bondlength(&self) -> f64 {
        self.target_bondlength
    }

    pub fn export_dir(&self) -> &PathBuf {
        &self.export_dir
    }

    pub fn kpoint_quality(&self) -> &KPointQuality {
        &self.kpoint_quality
    }

    pub fn edft(&self) -> bool {
        self.edft
    }

    pub fn potential_dir(&self) -> Option<&String> {
        self.potential_dir.as_ref()
    }
    pub fn x_range(&self) -> FractionalCoordRange {
        FractionalCoordRange::new(self.x_range.0, self.x_range.1)
    }
    pub fn y_range(&self) -> FractionalCoordRange {
        FractionalCoordRange::new(self.y_range.0, self.y_range.1)
    }
    pub fn z_range(&self) -> FractionalCoordRange {
        FractionalCoordRange::new(self.z_range.0, self.z_range.1)
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
