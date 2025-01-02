#![allow(dead_code)]
use std::fs::{self, create_dir};

use anyhow::Result;
use castep_cell_io::cell_document::CellDocument;
use chemrust_core::data::lattice::CrystalModel;
use clap::Parser;
use crystal_cif_io::DataBlock;

use castep_seeding::RootJobs;
use rhino_lib::{
    arg_parser::{Args, ProgramMode},
    interactive_ui::RunOptions,
    run, ModelFormat, RhinoExport, SearchJob, TaskTable,
};

fn main() -> Result<()> {
    let args = Args::parse();
    let program_mode = args.mode.unwrap_or(ProgramMode::I);
    match program_mode {
        ProgramMode::C => run_by_config(args.config_loc)?,
        ProgramMode::I => interactive_cli()?,
    }
    Ok(())
}

fn search_n_export<C: CrystalModel>(
    base_model: &C,
    search_job: &impl SearchJob,
    rhino_export: &impl RhinoExport,
) -> Result<()> {
    let results = run(base_model, search_job)?;
    rhino_export.export_all_kinds_sites::<CellDocument>(base_model, &results)?;
    rhino_export.export_all_kinds_sites::<DataBlock>(base_model, &results)?;
    Ok(())
}

fn run_by_table(yaml_table: &TaskTable) -> Result<()> {
    if !yaml_table.export_dir().exists() {
        create_dir(yaml_table.export_dir())?;
    }
    let model = ModelFormat::load_model(yaml_table.model_path())?;
    match &model {
        ModelFormat::Cell(cell_document) => {
            search_n_export(cell_document, yaml_table.search_config(), yaml_table)
        }
        ModelFormat::CifDataBlock(data_block) => {
            search_n_export(data_block, yaml_table.search_config(), yaml_table)
        }
    }?;
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
