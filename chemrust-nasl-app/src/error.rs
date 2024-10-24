use std::fmt::Display;

#[derive(Debug)]
pub enum RunError {
    Message(String),
    FormatError(FormatError),
    IO,
}

#[derive(Debug, Clone, Copy)]
pub enum FormatError {
    Identified,
    Supported,
    Compatible,
    ReadToString,
}

impl Display for FormatError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FormatError::Supported => f.write_str("This format is not supported"),
            FormatError::Compatible => {
                f.write_str("The file content is not compatible with the format requirement")
            }

            FormatError::Identified => {
                f.write_str("The file does not have extension suffix `.xxx`")
            }
            FormatError::ReadToString => f.write_str("Failed to read to string"),
        }
    }
}

impl Display for RunError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RunError::Message(m) => f.write_str(m),
            RunError::FormatError(m) => write!(f, "{m}"),
            RunError::IO => f.write_str("Error in IO operations"),
        }
    }
}

impl std::error::Error for RunError {}
impl std::error::Error for FormatError {}
