#![allow(dead_code, unused_imports)]
use std::{
    collections::HashMap,
    path::{Path, PathBuf},
};

use castep_cell_io::{CastepTask, CellDocument};
use castep_seeding::{CellBuilding, ParamBuilding};

use crate::{error::RunError, execution::ExportFile};
pub use potentials_io::{get_all_potentials, get_potential_entries, PotentialFileBytes};
pub use virt_files::{ContentStorage, FileModel};

mod potentials_io;

mod virt_files;

#[derive(Debug, Clone)]
pub struct Configurator {
    use_edft: bool,
    potentials_loc: PathBuf,
}

impl Configurator {
    pub fn new(use_edft: bool, potentials_loc: PathBuf) -> Self {
        Self {
            use_edft,
            potentials_loc,
        }
    }
}

impl CellBuilding for Configurator {}

impl ParamBuilding for Configurator {
    fn build_param_for_task(
        &self,
        template_cell: &castep_cell_io::CellDocument,
        castep_task: castep_cell_io::CastepTask,
    ) -> Result<castep_param_io::param::CastepParam, castep_seeding::SeedingErrors> {
        match castep_task {
            castep_cell_io::CastepTask::BandStructure => self.dos_param_template(
                template_cell,
                castep_cell_io::EnergyCutoff::Ultrafine,
                self.use_edft,
                &self.potentials_loc,
            ),
            castep_cell_io::CastepTask::GeometryOptimization => self.geom_opt_param_template(
                template_cell,
                castep_cell_io::EnergyCutoff::Ultrafine,
                self.use_edft,
                &self.potentials_loc,
            ),
        }
    }
}

pub fn create_seed_files<P: AsRef<Path>>(
    base_doc: &ExportFile<CellDocument, P>,
    config: &Configurator,
    potentials: &HashMap<String, PotentialFileBytes>,
) -> Result<ContentStorage, RunError> {
    let target_dir = Path::new(base_doc.file_name().as_ref().file_stem().unwrap());
    let potential_entries = get_potential_entries(base_doc)
        .iter()
        .map(|pot| {
            create_file(
                target_dir,
                pot.clone(),
                potentials.get(pot).unwrap().0.clone(),
            )
        })
        .collect::<Vec<FileModel>>();
    let file_stem = base_doc
        .file_name()
        .as_ref()
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap();
    let geom_cell = create_file(
        target_dir,
        base_doc.file_name(),
        config
            .build_cell_for_task(base_doc.file(), CastepTask::GeometryOptimization)
            .to_string()
            .into_bytes(),
    );
    let dos_cell = create_file(
        target_dir,
        format!("{}_DOS.cell", file_stem),
        config
            .build_cell_for_task(base_doc.file(), CastepTask::BandStructure)
            .to_string()
            .into_bytes(),
    );
    let geom_param = create_file(
        target_dir,
        format!("{}.param", file_stem),
        config
            .build_param_for_task(base_doc.file(), CastepTask::GeometryOptimization)?
            .to_string()
            .into_bytes(),
    );
    let dos_param = create_file(
        target_dir,
        format!("{}_DOS.param", file_stem),
        config
            .build_param_for_task(base_doc.file(), CastepTask::BandStructure)?
            .to_string()
            .into_bytes(),
    );
    Ok(ContentStorage::new(
        [
            potential_entries,
            [geom_cell, dos_cell, geom_param, dos_param].to_vec(),
        ]
        .concat(),
    ))
}

fn create_file<P: AsRef<Path>>(target_dir: &Path, filename: P, file: Vec<u8>) -> FileModel {
    FileModel::new(target_dir.join(filename), file)
}
