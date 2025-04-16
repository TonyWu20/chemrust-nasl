pub mod arg_parser;
pub mod error;
pub mod execution;
pub mod interactive_ui;
pub mod seeding;
pub mod supportive_data;
pub mod yaml_parser;

use chemrust_core::data::lattice::CrystalModel;
use chemrust_nasl::SearchReports;
use error::RunError;
pub use execution::{ExportFormat, ModelFormat, RhinoExport, SearchJob};
pub use interactive_ui::KPointQuality;
pub use yaml_parser::TaskTable;

pub fn run(
    model: &impl CrystalModel,
    search_job: &impl SearchJob,
) -> Result<SearchReports, RunError> {
    let results = search_job.search(model)?;
    let mul_empty = results.points().map(|v| v.is_empty()).unwrap_or(false);
    let single_empty = results
        .viable_single_points()
        .map(|v| v.is_empty())
        .unwrap_or(false);
    let double_empty = results
        .viable_double_points()
        .map(|v| v.is_empty())
        .unwrap_or(false);
    if mul_empty && single_empty && double_empty {
        return Err(RunError::NoAvailableResults);
    }
    Ok(results)
}
