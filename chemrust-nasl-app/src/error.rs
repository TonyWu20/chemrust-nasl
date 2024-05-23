use std::fmt::Display;

#[derive(Debug)]
pub enum RunError {
    Message(String),
}

impl Display for RunError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RunError::Message(m) => f.write_str(m),
        }
    }
}

impl std::error::Error for RunError {}
