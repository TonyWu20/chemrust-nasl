use std::{fs::read_to_string, path::Path};

use castep_cell_io::{CellDocument, CellParser};

use crate::error::FormatError;

pub fn load_cell_file<P: AsRef<Path>>(cell_path: P) -> Result<CellDocument, FormatError> {
    let content = read_to_string(cell_path).map_err(|_| FormatError::ReadToString)?;
    Ok(CellParser::from(&content)
        .parse()
        .map_err(|_| FormatError::Compatible)?)
}
