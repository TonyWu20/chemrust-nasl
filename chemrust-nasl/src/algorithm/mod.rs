use kiddo::ImmutableKdTree;
use nalgebra::Point3;

use self::{circle_check::check_circles, sphere_check::sphere_check};

use crate::coordination_sites::{CoordCircle, CoordPoint, CoordSphere};

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
pub struct SearchResults {
    spheres: Vec<CoordSphere>,
    circles: Vec<CoordCircle>,
    points: Vec<CoordPoint>,
}

impl SearchResults {
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
    SearchResults::new(
        sphere_intersect_results.spheres().to_vec(),
        pure_circles.to_vec(),
        dedup_points,
    )
}
