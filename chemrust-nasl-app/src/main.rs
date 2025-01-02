#![allow(dead_code)]
use std::fs::{self, create_dir};

use anyhow::Result;
use castep_cell_io::cell_document::CellDocument;
use clap::Parser;
use crystal_cif_io::DataBlock;

use arg_parser::{Args, ProgramMode};
use castep_seeding::RootJobs;
use interactive_ui::RunOptions;
use rhino_lib::{arg_parser, interactive_ui, run, ModelFormat, RhinoExport, TaskTable};

fn main() -> Result<()> {
    let args = Args::parse();
    let program_mode = args.mode.unwrap_or(ProgramMode::I);
    match program_mode {
        ProgramMode::C => run_by_config(args.config_loc)?,
        ProgramMode::I => interactive_cli()?,
    }
    Ok(())
}

fn run_by_table(yaml_table: &TaskTable) -> Result<()> {
    if !yaml_table.export_dir().exists() {
        create_dir(yaml_table.export_dir())?;
    }
    let model = ModelFormat::load_model(yaml_table.model_path())?;
    let results = match &model {
        ModelFormat::Cell(cell_document) => run(cell_document, yaml_table.search_config())?,
        ModelFormat::CifDataBlock(data_block) => run(data_block, yaml_table.search_config())?,
    };
    match model {
        ModelFormat::Cell(cell_document) => {
            yaml_table.export_all_kinds_sites::<CellDocument>(&cell_document, &results)?;
            yaml_table.export_all_kinds_sites::<DataBlock>(&cell_document, &results)?;
        }
        ModelFormat::CifDataBlock(data_block) => {
            yaml_table.export_all_kinds_sites::<CellDocument>(&data_block, &results)?;
            yaml_table.export_all_kinds_sites::<DataBlock>(&data_block, &results)?;
        }
    }
    if yaml_table.export_config().build_seed() {
        yaml_table.build_all(
            yaml_table.export_config(),
            yaml_table.export_config(),
            yaml_table.export_config().potential_loc(),
        )?;
        println!("Built all seed folders")
    }
    Ok(())
}

fn run_by_config(yaml_config_path: Option<String>) -> Result<()> {
    let filepath = yaml_config_path.unwrap_or("config.yaml".to_string());
    let yaml_table = TaskTable::load_task_table(filepath)?;
    run_by_table(&yaml_table)?;
    println!(
        "Results have been written to {}",
        yaml_table.export_dir().display()
    );
    Ok(())
}

fn interactive_cli() -> Result<()> {
    // CLI interpretation
    let run_options = RunOptions::new().unwrap();
    let yaml_table = run_options.build_task_table()?;
    run_by_table(&yaml_table)?;
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
