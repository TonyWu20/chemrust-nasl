use std::{fs::read_to_string, path::Path};

use chemrust_geometry::Structure;
use crystal_cif_io::parse_core_cif;

use crate::error::FormatError;
use crate::execution::model::cif_to_structure;

pub fn load_cif_file<P: AsRef<Path>>(cif_path: P) -> Result<Structure, FormatError> {
    let content = read_to_string(cif_path.as_ref()).map_err(|_| FormatError::ReadToString)?;
    let document = parse_core_cif(&content).map_err(|_| FormatError::Compatible)?;
    document
        .blocks
        .first()
        .cloned()
        .ok_or(FormatError::Missing)
        .and_then(|block| cif_to_structure(&block).map_err(|_| FormatError::Compatible))
}
