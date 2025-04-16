use std::io;

use castep_seeding::SeedingErrors;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum RunError {
    #[error("{0}")]
    Message(String),
    #[error("{0}")]
    FormatError(#[from] FormatError),
    #[error("Error in creating directory: {0} ")]
    CreateDir(#[from] io::Error),
    #[error("No available results. You may check if the atoms in the `.cell` are too close to the boundary of the lattice. Adjust them to be within the lattice could help.")]
    NoAvailableResults,
    #[error("Failed to load yaml config")]
    LoadYAML(#[from] serde_yaml::Error),
    #[error("Error in creating seed files: {0}")]
    SeedingError(#[from] SeedingErrors),
}

#[derive(Error, Debug, Clone, Copy)]
pub enum FormatError {
    #[error("The file does not have extension suffix `.xxx`")]
    Identified,
    #[error("This format is not supported")]
    Supported,
    #[error("The file content is not compatible with the designated format")]
    Compatible,
    #[error("Failed while reading to string")]
    ReadToString,
    #[error("Required parts are missing")]
    Missing,
}
