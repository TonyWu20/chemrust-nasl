use castep_periodic_table::element::ElementSymbol;
use serde::{Deserialize, Serialize};

use crate::supportive_data::FractionalCoordRange;

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct RunConfig {
    pub(crate) model_path: String,
    pub(crate) new_element: ElementSymbol,
    pub(crate) target_bondlength: f64,
    pub(crate) x_range: (f64, f64),
    pub(crate) y_range: (f64, f64),
    pub(crate) z_range: (f64, f64),
}

impl RunConfig {
    pub fn new(
        model_path: String,
        new_element: ElementSymbol,
        target_bondlength: f64,
        x_range: (f64, f64),
        y_range: (f64, f64),
        z_range: (f64, f64),
    ) -> Self {
        Self {
            model_path,
            new_element,
            target_bondlength,
            x_range,
            y_range,
            z_range,
        }
    }

    pub fn x_range(&self) -> FractionalCoordRange {
        FractionalCoordRange::new(self.x_range.0, self.x_range.1)
    }

    pub fn y_range(&self) -> FractionalCoordRange {
        FractionalCoordRange::new(self.y_range.0, self.y_range.1)
    }

    pub fn z_range(&self) -> FractionalCoordRange {
        FractionalCoordRange::new(self.z_range.0, self.z_range.1)
    }
}
