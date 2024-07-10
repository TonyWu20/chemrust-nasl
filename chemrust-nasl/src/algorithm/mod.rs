use std::{f64::consts::FRAC_PI_8, ops::ControlFlow};

use kiddo::SquaredEuclidean;
use nalgebra::{Point3, Vector3};
use rayon::prelude::*;

use self::sphere_check::SphereCheckResult;

use crate::{
    coordination_sites::{CoordCircle, MultiCoordPoint},
    geometry::{approx_cmp_f64, FloatOrdering},
    DelegatePoint,
};

mod circle_check;
mod search_config;
mod search_report;
mod site_index;
mod sphere_check;
#[cfg(test)]
mod test;

pub use search_config::SearchConfig;
pub use search_report::SearchReports;
pub use site_index::SiteIndex;

#[derive(Debug, Clone, Copy)]
pub struct AdsSiteLocator<'a> {
    site_index: &'a SiteIndex,
    config: &'a SearchConfig<'a>,
}

impl<'a> AdsSiteLocator<'a> {
    pub fn new(site_index: &'a SiteIndex, config: &'a SearchConfig<'a>) -> Self {
        Self { site_index, config }
    }
    pub fn search_sites(&self) -> SearchReports {
        let sphere_intersect_results = self.sphere_check();
        let special_sites = self.search_special_sites(&sphere_intersect_results);
        let viable_single_sites = self.search_possible_single_points();
        let viable_double_sites =
            self.search_possible_double_points(sphere_intersect_results.unchecked_circles());
        SearchReports::new(special_sites, viable_single_sites, viable_double_sites)
    }

    pub fn site_index(&self) -> &SiteIndex {
        self.site_index
    }

    pub fn config(&self) -> &SearchConfig<'a> {
        self.config
    }

    fn search_special_sites(
        &self,
        sphere_intersect_results: &SphereCheckResult,
    ) -> Option<Vec<MultiCoordPoint>> {
        let circle_check_results = self.check_circles(sphere_intersect_results.unchecked_circles());
        let points = [
            sphere_intersect_results.single_points(),
            circle_check_results.points(),
        ]
        .concat();
        let dedup_points = MultiCoordPoint::dedup_points(
            &points,
            self.site_index().coord_tree(),
            self.config().bondlength(),
        );
        if !dedup_points.is_empty() {
            println!("Special multi-coordinated sites search completed.");
            Some(dedup_points)
        } else {
            None
        }
    }
    fn search_possible_single_points(&self) -> Option<Vec<DelegatePoint<1>>> {
        let results: Vec<DelegatePoint<1>> = self
            .config()
            .to_check()
            .par_iter()
            .filter_map(|&(i, pt)| {
                self.brute_force_search_single(pt)
                    .map(|coord| DelegatePoint::<1>::new(coord, [i]))
            })
            .collect();
        if !results.is_empty() {
            Some(results)
        } else {
            None
        }
    }
    fn search_possible_double_points(
        &self,
        unchecked_circles: &[CoordCircle],
    ) -> Option<Vec<DelegatePoint<2>>> {
        let results: Vec<DelegatePoint<2>> = unchecked_circles
            .par_iter()
            .filter_map(|circ| {
                circ.get_possible_point(self.site_index().coord_tree(), self.config().bondlength())
            })
            .collect();
        if !results.is_empty() {
            Some(results)
        } else {
            None
        }
    }
    fn brute_force_search_single(&self, origin: Point3<f64>) -> Option<Point3<f64>> {
        let dist = self.config().bondlength();
        let kdtree = self.site_index().coord_tree();
        let step = FRAC_PI_8 / 2_f64;
        let azimuth: [f64; 32] = (0..32)
            .map(|i| i as f64 * step)
            .collect::<Vec<f64>>()
            .try_into()
            .unwrap();
        let elevation: [f64; 8] = (0..8)
            .map(|i| i as f64 * step)
            .collect::<Vec<f64>>()
            .try_into()
            .unwrap();
        let positions = elevation
            .iter()
            .map(|e| {
                azimuth
                    .iter()
                    .map(|a| {
                        let z = dist * e.sin();
                        let y = dist * e.cos() * a.sin();
                        let x = dist * e.cos() * a.cos();
                        Vector3::new(x, y, z)
                    })
                    .collect()
            })
            .collect::<Vec<Vec<Vector3<f64>>>>()
            .concat();
        let initial = Vector3::z_axis().scale(dist);
        let candidates = [[initial].to_vec(), positions].concat();
        let p = candidates.iter().try_for_each(|dir| {
            let p = origin + dir;
            if kdtree
                .within::<SquaredEuclidean>(&p.into(), (dist + 10_f64 * f64::EPSILON).powi(2))
                .iter()
                .any(|nb| {
                    matches!(
                        approx_cmp_f64(nb.distance, dist.powi(2)),
                        FloatOrdering::Less
                    )
                })
            {
                ControlFlow::Continue(())
            } else {
                ControlFlow::Break(p)
            }
        });
        match p {
            ControlFlow::Continue(_) => None,
            ControlFlow::Break(point) => Some(point),
        }
    }
}
