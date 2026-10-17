use castep_periodic_table::element::ElementSymbol;
use chemrust_core::data::lattice::{CrystalModel, UnitCellParameters};
use chemrust_nasl::{CoordSite, Visualize};
use crystal_cif_io::write_core_cif;
use std::path::Path;

use castep_cell_io::{
    cell_document::{CellDocument, IonicPosition},
    to_cell_document,
};

use crate::error::{FormatError, RunError};
use crate::execution::cif_model::{export_block, CifModel};

use super::{
    export::{CifExport, ExportFormat},
    format_loader::load_cif_file,
    load_cell_file,
};

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
    Cif(CifModel),
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

    pub fn as_cell(&self) -> Option<&CellDocument> {
        if let Self::Cell(v) = self {
            Some(v)
        } else {
            None
        }
    }

    pub fn as_cif(&self) -> Option<&CifModel> {
        if let Self::Cif(v) = self {
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
        element_symbol: &ElementSymbol,
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

impl ExportFormat for CifExport {
    type Item = String;
    fn add_new_site(
        base_model: &impl CrystalModel,
        coord_site: &(impl CoordSite + Visualize),
        element_symbol: &ElementSymbol,
    ) -> Self::Item {
        let new_name = format!("new_{}", coord_site.connecting_atoms_msg());
        let block = export_block(base_model, coord_site, element_symbol, &new_name);
        write_core_cif(std::slice::from_ref(&block))
    }

    fn suffix() -> String {
        "cif".to_string()
    }
}
