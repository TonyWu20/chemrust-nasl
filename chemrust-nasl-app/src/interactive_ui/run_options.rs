use std::path::Path;

use castep_periodic_table::element::ElementSymbol;
use inquire::{required, validator::Validation, CustomType, InquireError, Text};

use crate::{
    supportive_data::FractionalCoordRange,
    yaml_parser::{ExportConfig, RunConfig, TaskTable},
};

use super::{filepath_completer::FilePathCompleter, ExportOptions};

#[derive(Debug)]
pub struct RunOptions {
    filepath: String,
    new_element: ElementSymbol,
    target_bondlength: f64,
    x_range: FractionalCoordRange,
    y_range: FractionalCoordRange,
    z_range: FractionalCoordRange,
}

impl RunOptions {
    fn ask_filename() -> Result<String, InquireError> {
        let current_dir = std::env::current_dir().unwrap();
        let help_message = format!("Current directory: {}", current_dir.to_string_lossy());
        Text::new("Filepath of the model cell file:")
            .with_autocomplete(FilePathCompleter::default())
            .with_validator(required!("This field is required"))
            .with_validator(|input: &str| {
                if input.split('.').last().unwrap() == "cell" {
                    Ok(Validation::Valid)
                } else {
                    Ok(Validation::Invalid(
                        inquire::validator::ErrorMessage::Custom(
                            "Please enter the filepath of a `.cell` file".into(),
                        ),
                    ))
                }
            })
            .with_help_message(&help_message)
            .prompt()
    }
    fn ask_element() -> Result<ElementSymbol, InquireError> {
        CustomType::<ElementSymbol>::new("Element symbol of the new atom: ").prompt()
    }
    fn ask_bondlength() -> Result<f64, InquireError> {
        CustomType::<f64>::new("What is the target bondlength (Å)?").with_error_message("Please type a valid number").with_help_message("Type the desired bondlength between the new element atom and the existing atoms in model").prompt()
    }
    fn ask_frac_range(axis: &str) -> Result<FractionalCoordRange, InquireError> {
        let hint_message = format!(
            "Enter the fractional coordinate range to search in the direction of {}:",
            axis
        );
        let min = CustomType::<f64>::new(&hint_message)
            .with_help_message("Enter the lower limit, greater or equal to 0.0; press enter for default value (0.0)")
            .with_default(0.0)
            .prompt()?;
        let max = CustomType::<f64>::new(&hint_message)
            .with_help_message(
                "Enter the upper limit, less or equal to 1.0; press enter for default value (1.0)",
            )
            .with_default(1.0)
            .prompt()?;
        Ok(FractionalCoordRange::new(min, max))
    }
    pub fn new() -> Result<RunOptions, InquireError> {
        let filename = Self::ask_filename()?;
        let new_element = Self::ask_element()?;
        let target_bondlength = Self::ask_bondlength()?;
        let x_range = Self::ask_frac_range("x-axis")?;
        let y_range = Self::ask_frac_range("y-axis")?;
        let z_range = Self::ask_frac_range("z-axis")?;
        Ok(RunOptions {
            filepath: filename,
            new_element,
            target_bondlength,
            x_range,
            y_range,
            z_range,
        })
    }

    pub fn build_task_table(&self) -> Result<TaskTable, InquireError> {
        let run_config = RunConfig::new(
            self.filepath().into(),
            self.new_element,
            self.target_bondlength(),
            (self.x_range.min(), self.x_range.max()),
            (self.y_range.min(), self.y_range.max()),
            (self.z_range.min(), self.z_range.max()),
        );
        let export_config = self.export_config()?;
        Ok(TaskTable::new(run_config, export_config))
    }

    pub fn export_config(&self) -> Result<ExportConfig, InquireError> {
        let model_seedname = Path::new(self.filepath())
            .file_stem()
            .unwrap()
            .to_str()
            .unwrap();
        let export_options =
            ExportOptions::new(&self.new_element, self.target_bondlength(), model_seedname)?;
        let export_dir = export_options.export_dir().into();
        let potential_dir = Some(export_options.potential_dir().into());
        let kpoint_quality = export_options.kpoint_quality().to_owned();
        let edft = export_options.edft();
        let build_seed = export_options.build_seed();
        Ok(ExportConfig::new(
            export_dir,
            potential_dir,
            kpoint_quality,
            edft,
            build_seed,
        ))
    }

    pub fn filepath(&self) -> &str {
        self.filepath.as_ref()
    }

    pub fn new_element(&self) -> &ElementSymbol {
        &self.new_element
    }

    pub fn target_bondlength(&self) -> f64 {
        self.target_bondlength
    }
}
