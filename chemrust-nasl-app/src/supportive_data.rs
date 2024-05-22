use std::fmt::Display;

#[derive(Debug, Clone, Copy)]
pub struct FractionalCoordRange(f64, f64);

impl Display for FractionalCoordRange {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({}, {})", self.0, self.1)
    }
}

impl FractionalCoordRange {
    /// Creates a `FractionalCoordRange` with the range restricted to `[0.0, 1.0]`
    /// # Panics
    /// Panics if `lower > higher`
    pub fn new(lower: f64, higher: f64) -> Self {
        assert!(lower <= higher);
        Self(lower.clamp(0.0, 1.0), higher.clamp(0.0, 1.0))
    }
    pub fn is_in_range(&self, value: f64) -> bool {
        self.0 <= value && self.1 >= value
    }
    pub fn min(&self) -> f64 {
        self.0
    }
    pub fn max(&self) -> f64 {
        self.1
    }
}
