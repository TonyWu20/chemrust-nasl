use std::{error::Error, fmt::Display, path::Path};

#[derive(Debug, Clone, Copy)]
pub enum AcceptFormat {
    Cell,
}

#[derive(Debug, Clone, Copy)]
pub enum FormatError {
    NotIdentified,
    NotSupported,
    NotCompatible,
}

impl Display for FormatError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FormatError::NotSupported => f.write_str("This format is not supported"),
            FormatError::NotCompatible => {
                f.write_str("The file content is not compatible with the format requirement")
            }
            FormatError::NotIdentified => {
                f.write_str("The file does not have extension suffix `.xxx`")
            }
        }
    }
}

impl Error for FormatError {}

pub fn match_format<P: AsRef<Path>>(file_path: &P) -> Result<AcceptFormat, FormatError> {
    let suffix = file_path
        .as_ref()
        .extension()
        .ok_or(FormatError::NotIdentified)?
        .to_str()
        .unwrap();
    match suffix {
        "cell" => Ok(AcceptFormat::Cell),
        _ => Err(FormatError::NotSupported),
    }
}
