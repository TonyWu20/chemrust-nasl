use chemrust_geometry::Structure;
use std::path::Path;

use crate::error::{FormatError, RunError};

use super::{format_loader::load_cif_file, load_cell_file};

#[derive(Debug, Clone, Copy)]
pub enum AcceptFormat {
    Cell,
    Cif,
}

impl AcceptFormat {
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
}

/// A loaded model, as a `chemrust-geometry` `Structure`.
#[derive(Debug, Clone)]
pub enum ModelFormat {
    Cell(Structure),
    Cif(Structure),
}

impl ModelFormat {
    pub fn load_model<P: AsRef<Path>>(model_path: P) -> Result<ModelFormat, RunError> {
        let format = AcceptFormat::match_format(&model_path).map_err(RunError::FormatError)?;

        match format {
            AcceptFormat::Cell => Ok(ModelFormat::Cell(
                load_cell_file(model_path).map_err(RunError::FormatError)?,
            )),
            AcceptFormat::Cif => Ok(ModelFormat::Cif(
                load_cif_file(model_path).map_err(RunError::FormatError)?,
            )),
        }
    }

    pub fn as_model(&self) -> &Structure {
        match self {
            Self::Cell(v) | Self::Cif(v) => v,
        }
    }
}
