use std::{fs::read_to_string, path::Path};

use crystal_cif_io::{CifDocument, DataBlock};

use crate::error::FormatError;

pub fn load_cif_file<P: AsRef<Path>>(cif_path: P) -> Result<DataBlock, FormatError> {
    let content = read_to_string(cif_path.as_ref()).map_err(|_| FormatError::ReadToString)?;
    CifDocument::parse_from_str(&mut content.as_str())
        .map_err(|_| FormatError::Compatible)?
        .data_blocks()
        .ok_or(FormatError::Compatible)
        .and_then(|v| v.get(0).ok_or(FormatError::Missing))
        .cloned()
}
