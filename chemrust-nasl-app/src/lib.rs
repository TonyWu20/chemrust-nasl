pub mod arg_parser;
pub mod error;
pub mod execution;
pub mod interactive_ui;
pub mod supportive_data;
pub mod yaml_parser;

use error::RunError;
pub use interactive_ui::KPointQuality;
pub use yaml_parser::TaskTable;

pub fn run_by_table(task_table: &TaskTable) -> Result<(), RunError> {
    let results = execution::search(task_table)?;
    execution::export_results_in_cell(task_table, &results)?;
    println!(
        "Results have been written to {}",
        task_table.export_dir().display()
    );
    Ok(())
}
