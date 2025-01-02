use castep_periodic_table::element::ElementSymbol;
use chemrust_core::data::{
    atom::CoreAtomData,
    geom::coordinates::CoordData,
    lattice::{CrystalModel, UnitCellParameters},
};
use chemrust_nasl::{CoordSite, Visualize};
use std::path::Path;

use castep_cell_io::{
    cell_document::{CellDocument, IonicPosition},
    to_cell_document,
};
use crystal_cif_io::{
    data_dict::core_cif::{
        atom_site::chemrust_impl::from_atom_data, cell::chemrust_impl::from_unit_cell_parameters,
    },
    from_data_block_members, CifDocument, DataBlock,
};

use crate::error::{FormatError, RunError};

use super::{export::ExportFormat, format_loader::load_cif_file, load_cell_file};

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

#[derive(Debug, Clone)]
pub enum ModelFormat {
    Cell(CellDocument),
    CifDataBlock(DataBlock),
}

impl ModelFormat {
    pub fn load_model<P: AsRef<Path>>(model_path: P) -> Result<ModelFormat, RunError> {
        let format = AcceptFormat::match_format(&model_path).map_err(RunError::FormatError)?;

        match format {
            AcceptFormat::Cell => Ok(ModelFormat::Cell(
                load_cell_file(model_path).map_err(RunError::FormatError)?,
            )),
            AcceptFormat::Cif => Ok(ModelFormat::CifDataBlock(
                load_cif_file(model_path).map_err(RunError::FormatError)?,
            )),
        }
    }

    pub fn as_cell(&self) -> Option<&CellDocument> {
        if let Self::Cell(v) = self {
            Some(v)
        } else {
            None
        }
    }

    pub fn as_cif_data_block(&self) -> Option<&DataBlock> {
        if let Self::CifDataBlock(v) = self {
            Some(v)
        } else {
            None
        }
    }
}

impl ExportFormat for CellDocument {
    type Item = CellDocument;
    fn add_new_site(
        base_model: &impl CrystalModel,
        coord_site: &(impl CoordSite + Visualize),
        element_symbol: &castep_periodic_table::element::ElementSymbol,
    ) -> Self::Item {
        let mut new_cell: CellDocument = to_cell_document(base_model);
        let new_pos_coordinate =
            coord_site.fractional_coord(base_model.get_cell_parameters().lattice_bases());
        let new_pos = IonicPosition::new(*element_symbol, new_pos_coordinate.into(), None);
        new_cell
            .model_description_mut()
            .ionic_pos_block_mut()
            .positions_mut()
            .push(new_pos);
        new_cell
    }

    fn suffix() -> String {
        "cell".to_string()
    }
}

impl ExportFormat for DataBlock {
    type Item = CifDocument;
    fn add_new_site(
        base_model: &impl CrystalModel,
        coord_site: &(impl CoordSite + Visualize),
        element_symbol: &castep_periodic_table::element::ElementSymbol,
    ) -> Self::Item {
        let atom_data = base_model.get_atom_data();
        struct NewData {
            indices: Vec<usize>,
            symbols: Vec<ElementSymbol>,
            coords: Vec<CoordData>,
            labels: Vec<Option<String>>,
        }
        impl CoreAtomData for NewData {
            fn indices_repr(&self) -> Vec<usize> {
                self.indices.clone()
            }

            fn symbols_repr(&self) -> Vec<ElementSymbol> {
                self.symbols.clone()
            }

            fn coords_repr(&self) -> Vec<CoordData> {
                self.coords.clone()
            }

            fn labels_repr(&self) -> Vec<Option<String>> {
                self.labels.clone()
            }
        }
        let mut new_data = NewData {
            indices: atom_data.indices_repr(),
            symbols: atom_data.symbols_repr(),
            coords: atom_data.coords_repr(),
            labels: atom_data.labels_repr(),
        };
        let new_index = new_data.indices.len();
        new_data.indices.push(new_index);
        new_data.symbols.push(*element_symbol);
        new_data.coords.push(CoordData::Fractional(
            coord_site.fractional_coord(base_model.get_cell_parameters().lattice_bases()),
        ));
        new_data
            .labels
            .push(Some(format!("{}{}", element_symbol, new_index)));
        let new_name = format!("new_{}", coord_site.connecting_atoms_msg());
        let new_cell_block = from_unit_cell_parameters(base_model.get_cell_parameters());
        let new_data_block = from_atom_data(&new_data);
        let new_block = [new_cell_block, new_data_block].concat();
        from_data_block_members(&new_block, &new_name)
    }

    fn suffix() -> String {
        "cif".to_string()
    }
}
