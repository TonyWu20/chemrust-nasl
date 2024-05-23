use kiddo::{ImmutableKdTree, SquaredEuclidean};
use nalgebra::Point3;

use self::{circle_check::check_circles, sphere_check::sphere_check};

use crate::{
    coordination_sites::{CoordCircle, CoordPoint, CoordSphere},
    geometry::{approx_cmp_f64, FloatOrdering},
    Visualize,
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
    spheres: Vec<CoordSphere>,
    circles: Vec<CoordCircle>,
    points: Vec<CoordPoint>,
}

#[derive(Debug)]
pub enum SearchResults {
    Found(SearchReports),
    Empty,
}

impl SearchReports {
    pub fn new(
        spheres: Vec<CoordSphere>,
        circles: Vec<CoordCircle>,
        points: Vec<CoordPoint>,
    ) -> Self {
        Self {
            spheres,
            circles,
            points,
        }
    }

    pub fn spheres(&self) -> &[CoordSphere] {
        self.spheres.as_ref()
    }

    pub fn circles(&self) -> &[CoordCircle] {
        self.circles.as_ref()
    }

    pub fn points(&self) -> &[CoordPoint] {
        self.points.as_ref()
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

pub fn search_sites(site_index: &SiteIndex, search_config: &SearchConfig) -> SearchResults {
    let sphere_intersect_results = sphere_check(site_index, search_config);
    let circle_check_results = check_circles(
        sphere_intersect_results.unchecked_circles(),
        site_index,
        search_config,
    );
    let pure_circles = circle_check_results.circles();
    let points = [
        sphere_intersect_results.single_points(),
        circle_check_results.points(),
    ]
    .concat();
    let dedup_points =
        CoordPoint::dedup_points(&points, site_index.coord_tree(), search_config.bondlength);
    if sphere_intersect_results.spheres().is_empty()
        && pure_circles.is_empty()
        && dedup_points.is_empty()
    {
        SearchResults::Empty
    } else {
        SearchResults::Found(SearchReports::new(
            sphere_intersect_results.spheres().to_vec(),
            pure_circles.to_vec(),
            dedup_points,
        ))
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
