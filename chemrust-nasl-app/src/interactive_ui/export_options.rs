use std::{env, str::FromStr};

use castep_periodic_table::element::ElementSymbol;
use inquire::{Confirm, InquireError, Select, Text};

use super::{FilePathCompleter, KPointQuality, RunMode};

pub struct ExportOptions {
    export_dir: String,
    potential_dir: String,
    kpoint_quality: KPointQuality,
    edft: bool,
    build_seed: bool,
    run_mode: RunMode,
}

impl ExportOptions {
    pub fn new(
        new_element: &ElementSymbol,
        bondlength: f64,
        model_name: &str,
    ) -> Result<ExportOptions, InquireError> {
        let export_dir = Self::ask_export_dir(&new_element.to_string(), bondlength, model_name)?;
        let potential_dir = Self::ask_potential_dir()?;
        let kpoint_quality = Self::ask_kpoint_quality()?;
        let edft = Self::ask_edft(new_element)?;
        let build_seed = Self::ask_build_seed()?;
        let run_mode = Self::ask_run_mode()?;
        Ok(Self {
            export_dir,
            potential_dir,
            kpoint_quality,
            edft,
            build_seed,
            run_mode,
        })
    }
    fn ask_build_seed() -> Result<bool, InquireError> {
        Confirm::new("Do you want to generate castep job seed files?")
            .with_default(false)
            .prompt()
    }
    fn ask_export_dir(
        element_symbol: &str,
        bondlength: f64,
        model_name: &str,
    ) -> Result<String, InquireError> {
        Text::new("Please name the directory for exported seed files: ")
            .with_help_message("Default: element_bondlength_base-model-name")
            .with_autocomplete(FilePathCompleter::default())
            .with_default(&format!("{}_{}_{}", element_symbol, bondlength, model_name))
            .prompt()
    }
    fn ask_potential_dir() -> Result<String, InquireError> {
        let cwd = env::current_dir().expect("Current directory does not exist or no permissions");
        let potential_loc_path = format!("e.g.: {}/Potentials", cwd.display());
        Text::new("Please specify the location of castep psuedopotentials directory: ")
            .with_autocomplete(FilePathCompleter::default())
            .with_default(&potential_loc_path)
            .prompt()
    }
    fn ask_edft(new_element: &ElementSymbol) -> Result<bool, InquireError> {
        let edft_help_message = if (57..72).contains(&(*new_element as usize))
            || (89..104).contains(&(*new_element as usize))
        {
            format!(
                "The element {} belongs to the rare-earth series. edft method is suggested. (Type y/yes)",
                new_element
            )
        } else {
            format!(
                "dm method is suggested for the element {}. (Type n/no)",
                new_element
            )
        };
        Confirm::new("Use edft method for the `metals_method` option in CASTEP?(y/n or yes/no)")
            .with_help_message(&edft_help_message)
            .prompt()
    }
    fn ask_kpoint_quality() -> Result<KPointQuality, InquireError> {
        let kpoint_options: Vec<String> = vec!["Coarse".into(), "Medium".into(), "Fine".into()];
        let kpoint_quality =
            Select::new("Quality for k-point sampling?", kpoint_options).prompt()?;
        Ok(KPointQuality::from_str(&kpoint_quality).unwrap_or(KPointQuality::Coarse))
    }
    fn ask_run_mode() -> Result<RunMode, InquireError> {
        let run_mode_options: Vec<String> = vec![
            "Fast".into(),
            "Full".into(),
            "Dryrun".into(),
            "Debug".into(),
        ];
        let run_mode = Select::new("Run mode of program", run_mode_options).prompt()?;
        Ok(RunMode::from_str(&run_mode).unwrap_or(RunMode::Debug))
    }
    pub fn edft(&self) -> bool {
        self.edft
    }

    pub fn export_dir(&self) -> &str {
        self.export_dir.as_ref()
    }

    pub fn potential_dir(&self) -> &str {
        self.potential_dir.as_ref()
    }

    pub fn kpoint_quality(&self) -> &KPointQuality {
        &self.kpoint_quality
    }

    pub fn run_mode(&self) -> RunMode {
        self.run_mode
    }

    pub fn build_seed(&self) -> bool {
        self.build_seed
    }
}
