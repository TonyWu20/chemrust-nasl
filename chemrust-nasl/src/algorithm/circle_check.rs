use rayon::prelude::*;

use super::{SearchConfig, SiteIndex};

use crate::coordination_sites::{CoordCircle, CoordResult, MultiCoordPoint};

#[derive(Debug, Clone)]
pub struct CircleCheckResult {
    // circles: Vec<CoordCircle>,
    points: Vec<MultiCoordPoint>,
}

impl CircleCheckResult {
    pub fn new(points: Vec<MultiCoordPoint>) -> Self {
        Self { points }
    }

    // pub fn circles(&self) -> &[CoordCircle] {
    //     self.circles.as_ref()
    // }

    pub fn points(&self) -> &[MultiCoordPoint] {
        self.points.as_ref()
    }
}

pub fn check_circles(
    unchecked_circles: &[CoordCircle],
    site_index: &SiteIndex,
    search_config: &SearchConfig,
) -> CircleCheckResult {
    let kdtree = site_index.coord_tree();
    let points = site_index.coords();
    let dist = search_config.bondlength;
    // let mut coord_circles: Vec<CoordCircle> = Vec::new();
    let mut coord_points: Vec<MultiCoordPoint> = Vec::new();
    let check_results: Vec<CoordResult> = unchecked_circles
        .par_iter()
        .filter_map(|circ| -> Option<CoordResult> {
            circ.common_neighbours_intersect(kdtree, points, dist)
        })
        .collect();
    check_results.into_iter().for_each(|result| {
        if let CoordResult::Points(mut points) = result {
            coord_points.append(&mut points)
        }
    });
    CircleCheckResult::new(coord_points)
}
