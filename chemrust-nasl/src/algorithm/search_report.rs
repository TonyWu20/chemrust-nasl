use kiddo::SquaredEuclidean;

use crate::{
    approx_cmp_f64, AdsSiteLocator, DelegatePoint, FloatOrdering, MultiCoordPoint, Visualize,
};

#[derive(Debug)]
pub struct SearchReports {
    points: Option<Vec<MultiCoordPoint>>,
    viable_single_points: Option<Vec<DelegatePoint<1>>>,
    viable_double_points: Option<Vec<DelegatePoint<2>>>,
}

impl SearchReports {
    pub fn new(
        points: Option<Vec<MultiCoordPoint>>,
        viable_single_points: Option<Vec<DelegatePoint<1>>>,
        viable_double_points: Option<Vec<DelegatePoint<2>>>,
    ) -> Self {
        Self {
            points,
            viable_single_points,
            viable_double_points,
        }
    }

    pub fn points(&self) -> Option<&Vec<MultiCoordPoint>> {
        self.points.as_ref()
    }

    pub fn viable_single_points(&self) -> Option<&Vec<DelegatePoint<1>>> {
        self.viable_single_points.as_ref()
    }

    pub fn viable_double_points(&self) -> Option<&Vec<DelegatePoint<2>>> {
        self.viable_double_points.as_ref()
    }

    pub fn validated_results<T: Visualize + Clone>(
        coord_sites: &[T],
        site_locator: &AdsSiteLocator,
    ) -> Vec<T> {
        coord_sites
            .iter()
            .filter_map(|coord_site| validate_site(coord_site, site_locator))
            .cloned()
            .collect()
    }
}
pub fn validate_site<'a, T: Visualize>(
    coord_site: &'a T,
    site_locator: &AdsSiteLocator,
) -> Option<&'a T> {
    let coord = coord_site.determine_coord();
    let bondlength = site_locator.config().bondlength();
    let dist = bondlength.powi(2);
    if site_locator
        .site_index()
        .coord_tree()
        .within::<SquaredEuclidean>(&coord.into(), dist)
        .iter()
        .any(|nb| matches!(approx_cmp_f64(nb.distance, dist), FloatOrdering::Less))
    {
        None
    } else {
        Some(coord_site)
    }
}
