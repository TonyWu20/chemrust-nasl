#![allow(dead_code)]
use crate::coordination_sites::{CoordCircle, MultiCoordPoint};

#[derive(Debug, Clone)]
pub struct CircleCheckResult {
    circles: Vec<CoordCircle>,
    points: Vec<MultiCoordPoint>,
}

impl CircleCheckResult {
    pub fn new(circles: Vec<CoordCircle>, points: Vec<MultiCoordPoint>) -> Self {
        Self { circles, points }
    }

    pub fn circles(&self) -> &[CoordCircle] {
        self.circles.as_ref()
    }

    pub fn points(&self) -> &[MultiCoordPoint] {
        self.points.as_ref()
    }
}
