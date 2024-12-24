use std::path::Path;

use castep_cell_io::cell_document::CellDocument;
use crystal_cif_io::DataBlock;

use crate::error::FormatError;

#[derive(Debug, Clone, Copy)]
pub enum AcceptFormat {
    Cell,
    Cif,
}

pub enum ModelFormat {
    Cell(CellDocument),
    CifDataBlock(DataBlock),
}

pub fn match_format<P: AsRef<Path>>(file_path: P) -> Result<AcceptFormat, FormatError> {
    let suffix = file_path
        .as_ref()
        .extension()
        .ok_or(FormatError::Identified)?
        .to_str()
        .unwrap();
    match suffix {
        "cell" => Ok(AcceptFormat::Cell),
        "cif" => Ok(AcceptFormat::Cif),
        _ => Err(FormatError::Supported),
    }
}
