use std::{
    collections::{HashMap, HashSet},
    fs::read_to_string,
    path::Path,
};

use castep_cell_io::{
    cell_document::{
        sections::species_characters::{SpeciesBlock, SpeciesEntry},
        CellEntries,
    },
    CellDocument,
};
use castep_periodic_table::{
    data::ELEMENT_TABLE,
    element::{ElementSymbol, LookupElement},
};

use crate::{
    error::{FormatError, RunError},
    execution::ExportFile,
};

#[derive(Debug, Clone)]
pub struct PotentialFileBytes(pub(crate) Vec<u8>);

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

fn get_all_elements<P: AsRef<Path>>(
    base_cells: &[ExportFile<CellDocument, P>],
) -> HashSet<ElementSymbol> {
    HashSet::from_iter(
        base_cells
            .iter()
            .flat_map(|cell_doc| -> Vec<ElementSymbol> { cell_doc.file().get_elements() })
            .collect::<Vec<ElementSymbol>>(),
    )
}

pub fn get_potential_entries<P: AsRef<Path>>(
    cell_file: &ExportFile<CellDocument, P>,
) -> Vec<String> {
    cell_file
        .file()
        .other_entries()
        .and_then(|v| {
            v.iter()
                .find(|entry| matches!(entry, CellEntries::SpeciesPot(_sp)))
                .and_then(|entry| {
                    if let CellEntries::SpeciesPot(sp) = entry {
                        Some(
                            sp.items()
                                .iter()
                                .map(|s| s.item())
                                .cloned()
                                .collect::<Vec<String>>(),
                        )
                    } else {
                        None
                    }
                })
        })
        .unwrap_or(
            cell_file
                .file()
                .get_elements()
                .iter()
                .map(|elm| ELEMENT_TABLE.get_by_symbol(*elm).potential().into())
                .collect::<Vec<String>>(),
        )
}
