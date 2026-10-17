use std::{fs::read_to_string, path::Path};

use crystal_cif_io::parse_core_cif;

use crate::{error::FormatError, execution::cif_model::CifModel};

pub fn load_cif_file<P: AsRef<Path>>(cif_path: P) -> Result<CifModel, FormatError> {
    let content = read_to_string(cif_path.as_ref()).map_err(|_| FormatError::ReadToString)?;
    let document = parse_core_cif(&content).map_err(|_| FormatError::Compatible)?;
    document
        .blocks
        .first()
        .cloned()
        .ok_or(FormatError::Missing)
        .and_then(|block| CifModel::from_block(&block).map_err(|_| FormatError::Compatible))
}
