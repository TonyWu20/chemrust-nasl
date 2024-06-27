use std::{error::Error, fmt::Display, path::Path};

use castep_cell_io::CellDocument;
use chemrust_core::data::lattice::CrystalModel;

#[derive(Debug, Clone, Copy)]
pub enum AcceptFormat {
    Cell,
}

pub enum ModelFormat {
    Cell(CellDocument),
}

pub struct Model<T: CrystalModel>(pub(crate) T);

#[derive(Debug, Clone, Copy)]
pub enum FormatError {
    Identified,
    Supported,
    Compatible,
}

impl Display for FormatError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FormatError::Supported => f.write_str("This format is not supported"),
            FormatError::Compatible => {
                f.write_str("The file content is not compatible with the format requirement")
            }
            FormatError::Identified => {
                f.write_str("The file does not have extension suffix `.xxx`")
            }
        }
    }
}

impl Error for FormatError {}

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
