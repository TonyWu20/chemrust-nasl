#![allow(dead_code)]
use std::{error::Error, fs};

use arg_parser::ProgramMode;
use clap::Parser;
use error::RunError;
use interactive_ui::RunOptions;

use crate::{
    execution::{export_results_in_cell, search},
    yaml_parser::TaskTable,
};

mod arg_parser;
mod error;
mod execution;
mod interactive_ui;
mod supportive_data;
mod yaml_parser;

fn main() -> Result<(), Box<dyn Error>> {
    let args = arg_parser::Args::parse();
    let program_mode = args.mode.unwrap_or(ProgramMode::I);
    if program_mode == ProgramMode::I {
        interactive_cli()?;
    } else {
        run_by_config(args.config_loc)?;
    }
    Ok(())
}

fn run_by_config(yaml_config_path: Option<String>) -> Result<(), Box<dyn Error>> {
    let filepath = yaml_config_path.unwrap_or("config.yaml".to_string());
    let yaml_table = TaskTable::load_task_table(filepath)?;
    let results = search(&yaml_table)?;
    match results {
        chemrust_nasl::SearchResults::Found(report) => {
            println!(
                "Successful. Results have been written to {}",
                yaml_table.export_dir()
            );
            export_results_in_cell(&yaml_table, &report)
        }
        chemrust_nasl::SearchResults::Empty => {
            Err(RunError::Message("No results for this config".to_string()))?
        }
    }
}

fn interactive_cli() -> Result<(), Box<dyn Error>> {
    // CLI interpretation
    let run_options = RunOptions::new().unwrap();
    let yaml_table = run_options.export_config()?;
    let results = search(&yaml_table)?;
    match results {
        chemrust_nasl::SearchResults::Found(report) => {
            export_results_in_cell(&yaml_table, &report)?
        }
        chemrust_nasl::SearchResults::Empty => {
            Err(RunError::Message("No results for this config".to_string()))?
        }
    }
    let export_table_filename = format!(
        "{}/{}.yaml",
        yaml_table.export_dir(),
        yaml_table.export_dir()
    );
    fs::write(export_table_filename, serde_yaml::to_string(&yaml_table)?)?;
    Ok(())
}
