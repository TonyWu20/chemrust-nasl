#![allow(dead_code)]
use std::{error::Error, fs};

use clap::Parser;
use rhino_lib::arg_parser::Args;
use rhino_lib::arg_parser::ProgramMode;
use rhino_lib::interactive_ui::RunOptions;

use rhino_lib::execution::{export_results_in_cell, search};
use rhino_lib::run_by_table;
use rhino_lib::yaml_parser::TaskTable;

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();
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
    run_by_table(&yaml_table)?;
    Ok(())
}

fn interactive_cli() -> Result<(), Box<dyn Error>> {
    // CLI interpretation
    let run_options = RunOptions::new().unwrap();
    let yaml_table = run_options.export_config()?;
    let results = search(&yaml_table)?;
    let (mul, sing, doub) = export_results_in_cell(&yaml_table, &results)?;
    let export_table_filename = yaml_table.export_dir().join(
        yaml_table
            .export_dir()
            .file_name()
            .expect("ends with '..'")
            .to_str()
            .expect("Invalid Unicode"),
    );
    if mul == 0 && sing == 0 && doub == 0 {
        println!("No avaliable results. You may check if the atoms in the `.cell` are too close to the boundary of the lattice. Adjust them to be within the lattice could help.");
    } else {
        println!(
            "Results have been written to {}",
            yaml_table.export_dir().display()
        );
        fs::write(
            format!("{}.yaml", export_table_filename.display()),
            serde_yaml::to_string(&yaml_table)?,
        )?;
    }
    Ok(())
}
