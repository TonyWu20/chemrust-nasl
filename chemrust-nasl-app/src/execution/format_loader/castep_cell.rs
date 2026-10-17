use std::{fs::read_to_string, path::Path};

use castep_cell_fmt::parse;
use castep_cell_io::CellDocument;
use chemrust_geometry::Structure;

use crate::error::FormatError;
use crate::execution::model::cell_to_structure;

pub fn load_cell_file<P: AsRef<Path>>(cell_path: P) -> Result<Structure, FormatError> {
    let content = read_to_string(cell_path).map_err(|_| FormatError::ReadToString)?;
    load_cell_content(content)
}

pub fn load_cell_content(cell_content: String) -> Result<Structure, FormatError> {
    let doc: CellDocument = parse(&cell_content).map_err(|_| FormatError::Compatible)?;
    cell_to_structure(&doc).map_err(|_| FormatError::Compatible)
}
