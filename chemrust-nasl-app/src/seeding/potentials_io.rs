use std::{
    collections::{HashMap, HashSet},
    fs::read_to_string,
    path::Path,
};

use castep_cell_io::CellDocument;
use castep_periodic_table::{data::ELEMENT_TABLE, element::LookupElement};
use castep_seeding::cell_elements;

use crate::{
    error::{FormatError, RunError},
    execution::ExportFile,
};

#[derive(Debug, Clone)]
pub struct PotentialFileBytes(pub Vec<u8>);

pub fn get_all_potentials<P: AsRef<Path>, Q: AsRef<Path>>(
    base_cells: &[ExportFile<CellDocument, P>],
    potential_loc: Q,
) -> Result<HashMap<String, PotentialFileBytes>, RunError> {
    Ok(HashMap::<String, PotentialFileBytes>::from_iter(
        HashSet::<String>::from_iter(base_cells.iter().flat_map(get_potential_entries))
            .iter()
            .map(|potential_name| {
                Ok((
                    potential_name.clone(),
                    PotentialFileBytes(
                        read_to_string(potential_loc.as_ref().join(potential_name))
                            .map_err(|_| RunError::FormatError(FormatError::ReadToString))?
                            .into_bytes(),
                    ),
                ))
            })
            .collect::<Result<Vec<(String, PotentialFileBytes)>, RunError>>()?,
    ))
}

pub fn get_potential_entries<P: AsRef<Path>>(
    cell_file: &ExportFile<CellDocument, P>,
) -> HashSet<String> {
    let doc = cell_file.file();
    let species_pots = doc.species.species_pot.as_ref().map(|sp| {
        sp.potentials
            .iter()
            .map(|entry| entry.filename.clone())
            .collect::<Vec<String>>()
    });
    HashSet::from_iter(species_pots.unwrap_or_else(|| {
        cell_elements(doc)
            .iter()
            .map(|elm| ELEMENT_TABLE.get_by_symbol(*elm).potential().into())
            .collect::<Vec<String>>()
    }))
}
