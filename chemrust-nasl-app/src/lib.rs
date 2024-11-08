pub mod arg_parser;
pub mod error;
pub mod execution;
pub mod interactive_ui;
pub mod supportive_data;
pub mod yaml_parser;

use chemrust_nasl::SearchReports;
use error::RunError;
pub use interactive_ui::KPointQuality;
pub use yaml_parser::TaskTable;

pub fn report(results: &SearchReports) {
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

pub fn run_by_table(task_table: &TaskTable) -> Result<(), RunError> {
    let results = execution::search(task_table)?;
    report(&results);
    execution::export_results_in_cell(task_table, &results)?;
    println!(
        "Results have been written to {}",
        task_table.export_dir().display()
    );
    Ok(())
}
