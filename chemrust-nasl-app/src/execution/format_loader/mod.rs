use std::{error::Error, fs::read_to_string, path::Path};

use castep_cell_io::{CellDocument, CellParser};

pub fn load_cell_file<P: AsRef<Path>>(cell_path: P) -> Result<CellDocument, Box<dyn Error>> {
    let content = read_to_string(cell_path)?;
    Ok(CellParser::from(&content).parse()?)
}
