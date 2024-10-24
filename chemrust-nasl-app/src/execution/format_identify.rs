use std::path::Path;

use castep_cell_io::CellDocument;
use chemrust_core::data::lattice::CrystalModel;

use crate::error::FormatError;

#[derive(Debug, Clone, Copy)]
pub enum AcceptFormat {
    Cell,
}

pub enum ModelFormat {
    Cell(CellDocument),
}

pub struct Model<T: CrystalModel>(pub(crate) T);

pub fn match_format<P: AsRef<Path>>(file_path: &P) -> Result<AcceptFormat, FormatError> {
    let suffix = file_path
        .as_ref()
        .extension()
        .ok_or(FormatError::Identified)?
        .to_str()
        .unwrap();
    match suffix {
        "cell" => Ok(AcceptFormat::Cell),
        _ => Err(FormatError::Supported),
    }
}
