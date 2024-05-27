use std::{
    f64::{consts::FRAC_PI_8, EPSILON},
    ops::ControlFlow,
};

use kiddo::{ImmutableKdTree, SquaredEuclidean};
use nalgebra::{Point3, Vector3};
use rayon::prelude::*;

use self::{
    circle_check::check_circles,
    sphere_check::{sphere_check, SphereCheckResult},
};

use crate::{
    coordination_sites::{CoordCircle, MultiCoordPoint},
    geometry::{approx_cmp_f64, FloatOrdering},
    DelegatePoint, Visualize,
};

mod circle_check;
mod sphere_check;

#[derive(Debug, Clone, Copy)]
pub struct SearchConfig<'a> {
    to_check: &'a [(usize, Point3<f64>)],
    bondlength: f64,
}

impl<'a> SearchConfig<'a> {
    pub fn new(to_check: &'a [(usize, Point3<f64>)], bondlength: f64) -> Self {
        Self {
            to_check,
            bondlength,
        }
    }

    pub fn to_check(&self) -> &[(usize, Point3<f64>)] {
        self.to_check
    }

    pub fn bondlength(&self) -> f64 {
        self.bondlength
    }
}

#[cfg(test)]
mod test;

fn build_kd_tree_from_points(points: &[Point3<f64>]) -> ImmutableKdTree<f64, 3> {
    let entries = points.iter().map(|&p| p.into()).collect::<Vec<[f64; 3]>>();
    ImmutableKdTree::new_from_slice(&entries)
}

pub struct SiteIndex {
    coords: Vec<Point3<f64>>,
    coord_tree: ImmutableKdTree<f64, 3>,
}

impl SiteIndex {
    pub fn new(coords: Vec<Point3<f64>>) -> Self {
        let coord_tree = build_kd_tree_from_points(&coords);
        Self { coords, coord_tree }
    }

    pub fn coords(&self) -> &[Point3<f64>] {
        self.coords.as_ref()
    }

    pub fn coord_tree(&self) -> &ImmutableKdTree<f64, 3> {
        &self.coord_tree
    }
}

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
        site_index: &SiteIndex,
        search_config: &SearchConfig,
    ) -> Vec<T> {
        coord_sites
            .iter()
            .filter_map(|coord_site| validate_site(coord_site, site_index, search_config))
            .cloned()
            .collect()
    }
}

pub fn search_sites(site_index: &SiteIndex, search_config: &SearchConfig) -> SearchReports {
    let sphere_intersect_results = sphere_check(site_index, search_config);
    let special_sites = search_special_sites(&sphere_intersect_results, site_index, search_config);
    let viable_single_sites = search_possible_single_points(site_index, search_config);
    let viable_double_sites = search_possible_double_points(
        sphere_intersect_results.unchecked_circles(),
        site_index.coord_tree(),
        search_config.bondlength(),
    );
    SearchReports::new(special_sites, viable_single_sites, viable_double_sites)
}

fn search_special_sites(
    sphere_intersect_results: &SphereCheckResult,
    site_index: &SiteIndex,
    search_config: &SearchConfig,
) -> Option<Vec<MultiCoordPoint>> {
    let circle_check_results = check_circles(
        sphere_intersect_results.unchecked_circles(),
        site_index,
        search_config,
    );
    let points = [
        sphere_intersect_results.single_points(),
        circle_check_results.points(),
    ]
    .concat();
    let dedup_points =
        MultiCoordPoint::dedup_points(&points, site_index.coord_tree(), search_config.bondlength);
    if !dedup_points.is_empty() {
        println!("Special multi-coordinated sites search completed.");
        Some(dedup_points)
    } else {
        None
    }
}

fn search_possible_single_points(
    site_index: &SiteIndex,
    search_config: &SearchConfig,
) -> Option<Vec<DelegatePoint<1>>> {
    let results: Vec<DelegatePoint<1>> = search_config
        .to_check()
        .par_iter()
        .filter_map(|&(i, pt)| {
            brute_force(pt, search_config.bondlength(), site_index.coord_tree())
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
    unchecked_circles: &[CoordCircle],
    kdtree: &ImmutableKdTree<f64, 3>,
    dist: f64,
) -> Option<Vec<DelegatePoint<2>>> {
    let results: Vec<DelegatePoint<2>> = unchecked_circles
        .par_iter()
        .filter_map(|circ| circ.get_possible_point(kdtree, dist))
        .collect();
    if !results.is_empty() {
        Some(results)
    } else {
        None
    }
}

fn brute_force(
    origin: Point3<f64>,
    dist: f64,
    kdtree: &ImmutableKdTree<f64, 3>,
) -> Option<Point3<f64>> {
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
            .within::<SquaredEuclidean>(&p.into(), (dist + 10_f64 * EPSILON).powi(2))
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

pub fn validate_site<'a, 'b, T: Visualize>(
    coord_site: &'a T,
    site_index: &'b SiteIndex,
    search_config: &'b SearchConfig,
) -> Option<&'a T> {
    let coord = coord_site.determine_coord();
    let bondlength = search_config.bondlength;
    let dist = bondlength.powi(2);
    if site_index
        .coord_tree
        .within::<SquaredEuclidean>(&coord.into(), dist)
        .iter()
        .any(|nb| matches!(approx_cmp_f64(nb.distance, dist), FloatOrdering::Less))
    {
        None
    } else {
        Some(coord_site)
    }
}
