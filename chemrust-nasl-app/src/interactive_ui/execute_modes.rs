use std::{fmt::Display, str::FromStr};

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum RunMode {
    /// Generate seed files without copying potentials
    Fast,
    /// Generate seed files and copy potentials
    Full,
    /// Copy potentials after seed files generation
    Post,
    /// Reorganize
    Dryrun,
    /// Debug
    Debug,
    /// Clean the generated folder
    Clean,
}

impl Display for RunMode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RunMode::Fast => write!(f, "Fast"),
            RunMode::Full => write!(f, "Full"),
            RunMode::Post => write!(f, "Post"),
            RunMode::Dryrun => write!(f, "Dryrun"),
            RunMode::Debug => write!(f, "Debug"),
            RunMode::Clean => write!(f, "Clean"),
        }
    }
}

#[derive(Debug, PartialEq, Eq)]
pub struct ParseRunModeError;

impl FromStr for RunMode {
    type Err = ParseRunModeError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "Fast" => Ok(RunMode::Fast),
            "Full" => Ok(RunMode::Full),
            "Post" => Ok(RunMode::Post),
            "Dryrun" => Ok(RunMode::Dryrun),
            "Debug" => Ok(RunMode::Debug),
            "Clean" => Ok(RunMode::Clean),
            _ => Err(ParseRunModeError),
        }
    }
}
