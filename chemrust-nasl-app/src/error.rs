use thiserror::Error;

#[derive(Error, Debug)]
pub enum RunError {
    #[error("{0}")]
    Message(String),
    #[error("{0}")]
    FormatError(#[from] FormatError),
    #[error("Error in IO operations")]
    IO,
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
