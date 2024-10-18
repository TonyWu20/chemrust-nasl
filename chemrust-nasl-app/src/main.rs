#![allow(dead_code)]
use std::{error::Error, fs};

use arg_parser::ProgramMode;
use chemrust_nasl::SearchReports;
use clap::Parser;
use interactive_ui::RunOptions;

use crate::execution::{export_results_in_cell, search};

pub use yaml_parser::TaskTable;

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

fn report(results: &SearchReports) {
    println!(
        "Found {} multi-coordinated positions;",
        results.points().map(|v| v.len()).unwrap_or_default()
    );
    println!(
        "Found {} possible doubly-coordinated positions;",
        results.viable_double_points().map(|v| v.len()).unwrap_or(0)
    );
    println!(
        "Found {} possible singly-coordinated positions;",
        results.viable_single_points().map(|v| v.len()).unwrap_or(0)
    );
}

pub fn run_by_table(task_table: &TaskTable) -> Result<(), Box<dyn Error>> {
    let results = search(&task_table)?;
    report(&results);
    export_results_in_cell(&task_table, &results)?;
    println!(
        "Results have been written to {}",
        task_table.export_dir().display()
    );
    Ok(())
}

fn run_by_config(yaml_config_path: Option<String>) -> Result<(), Box<dyn Error>> {
    let filepath = yaml_config_path.unwrap_or("config.yaml".to_string());
    let yaml_table = TaskTable::load_task_table(filepath)?;
    run_by_table(&yaml_table)
}

fn interactive_cli() -> Result<(), Box<dyn Error>> {
    // CLI interpretation
    let run_options = RunOptions::new().unwrap();
    let yaml_table = run_options.export_config()?;
    let results = search(&yaml_table)?;
    report(&results);
    export_results_in_cell(&yaml_table, &results)?;
    let export_table_filename = yaml_table.export_dir().join(
        yaml_table
            .export_dir()
            .file_name()
            .expect("ends with '..'")
            .to_str()
            .expect("Invalid Unicode"),
    );
    println!(
        "Results have been written to {}",
        yaml_table.export_dir().display()
    );
    fs::write(
        format!("{}.yaml", export_table_filename.display()),
        serde_yaml::to_string(&yaml_table)?,
    )?;
    Ok(())
}
