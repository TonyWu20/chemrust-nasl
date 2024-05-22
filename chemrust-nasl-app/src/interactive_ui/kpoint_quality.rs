use std::{fmt::Display, str::FromStr};

use serde::{Deserialize, Serialize};

#[derive(Debug, Serialize, Deserialize, Clone)]
pub enum KPointQuality {
    Coarse,
    Medium,
    Fine,
}

impl Display for KPointQuality {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            KPointQuality::Coarse => write!(f, "Coarse"),
            KPointQuality::Medium => write!(f, "Medium"),
            KPointQuality::Fine => write!(f, "Fine"),
        }
    }
}

#[derive(Debug, PartialEq, Eq)]
pub struct ParseKPointQualityError;

impl FromStr for KPointQuality {
    type Err = ParseKPointQualityError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "Coarse" => Ok(KPointQuality::Coarse),
            "Medium" => Ok(KPointQuality::Medium),
            "Fine" => Ok(KPointQuality::Fine),
            _ => Err(ParseKPointQualityError),
        }
    }
}
