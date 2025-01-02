use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use crate::KPointQuality;

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct ExportConfig {
    pub(crate) export_dir: PathBuf,
    pub(crate) potential_dir: Option<PathBuf>,
    pub(crate) kpoint_quality: KPointQuality,
    pub(crate) edft: bool,
    #[serde(default)]
    build_seed: bool,
}

impl ExportConfig {
    pub fn new(
        export_dir: PathBuf,
        potential_dir: Option<PathBuf>,
        kpoint_quality: KPointQuality,
        edft: bool,
        build_seed: bool,
    ) -> Self {
        Self {
            export_dir,
            potential_dir,
            kpoint_quality,
            edft,
            build_seed,
        }
    }

    pub fn potential_loc(&self) -> &Path {
        self.potential_dir
            .as_ref()
            .map_or(Path::new("Potentials"), |v| v)
    }

    pub fn build_seed(&self) -> bool {
        self.build_seed
    }
}
