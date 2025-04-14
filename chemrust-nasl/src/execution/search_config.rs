use std::{f64::consts::FRAC_PI_8, ops::ControlFlow};

use kd_tree::KdIndexTree;
use nalgebra::{distance_squared, Point3, Vector3};
use rayon::prelude::*;

use crate::algorithm::{
    circle_check::CircleCheckResult,
    sphere_check::{sphere_check_fn, SphereCheckResult},
};

use crate::{
    coordination_sites::{CoordCircle, MultiCoordPoint},
    geometry::{approx_cmp_f64, FloatOrdering},
    CoordResult, DelegatePoint, Visualize,
};

#[derive(Debug, Clone)]
pub struct SearchConfig<'a> {
    to_check: &'a [(usize, Point3<f64>)],
    coord_tree: KdIndexTree<'a, Point3<f64>>,
    bondlength: f64,
}

impl<'a> SearchConfig<'a> {
    pub fn new(
        to_check: &'a [(usize, Point3<f64>)],
        all_points: &'a [Point3<f64>],
        bondlength: f64,
    ) -> Self {
        Self {
            to_check,
            coord_tree: KdIndexTree::build_by_ordered_float(all_points),
            bondlength,
        }
    }

    pub fn to_check(&self) -> &[(usize, Point3<f64>)] {
        self.to_check
    }

    pub fn bondlength(&self) -> f64 {
        self.bondlength
    }

    pub fn coord_tree(&self) -> &KdIndexTree<'a, Point3<f64>> {
        &self.coord_tree
    }

    fn check_circles(&self, unchecked_circles: &[CoordCircle]) -> CircleCheckResult {
        let kdtree = self.coord_tree();
        let points = self.coord_tree().source();
        let dist = self.bondlength;
        let mut coord_circles: Vec<CoordCircle> = Vec::new();
        let mut coord_points: Vec<MultiCoordPoint> = Vec::new();
        let check_results: Vec<CoordResult> = unchecked_circles
            .par_iter()
            .filter_map(|circ| -> Option<CoordResult> {
                circ.common_neighbours_intersect(kdtree, points, dist)
            })
            .collect();
        check_results.into_iter().for_each(|result| match result {
            CoordResult::Circle(c) => coord_circles.push(c),
            CoordResult::Points(mut points) => coord_points.append(&mut points),
            _ => (),
        });
        CircleCheckResult::new(coord_circles, coord_points)
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
        let dedup_points =
            MultiCoordPoint::dedup_points(&points, self.coord_tree(), self.bondlength);
        if !dedup_points.is_empty() {
            println!("Special multi-coordinated sites search completed.");
            Some(dedup_points)
        } else {
            None
        }
    }

    fn search_possible_single_points(&self) -> Option<Vec<DelegatePoint<1>>> {
        let results: Vec<DelegatePoint<1>> = self
            .to_check()
            .par_iter()
            .filter_map(|&(i, pt)| {
                brute_force(pt, self.bondlength(), self.coord_tree())
                    .map(|coord| DelegatePoint::<1>::new(coord, [i]))
            })
            .collect();
        if !results.is_empty() {
            Some(results)
        } else {
            None
        }
    }

    /// Conduct the search
    /// Validate the sites before output
    pub fn search_sites(&self) -> SearchReports {
        let sphere_intersect_results = self.sphere_check();
        let special_sites = self
            .search_special_sites(&sphere_intersect_results)
            .map(|v| {
                v.into_iter()
                    .filter_map(|site| self.validate_site(site))
                    .collect()
            });
        let viable_single_sites = self.search_possible_single_points().map(|v| {
            v.into_iter()
                .filter_map(|site| self.validate_site(site))
                .collect()
        });
        let viable_double_sites = self
            .search_possible_double_points(sphere_intersect_results.unchecked_circles())
            .map(|v| {
                v.into_iter()
                    .filter_map(|site| self.validate_site(site))
                    .collect()
            });
        SearchReports::new(special_sites, viable_single_sites, viable_double_sites)
    }

    fn search_possible_double_points(
        &self,
        unchecked_circles: &[CoordCircle],
    ) -> Option<Vec<DelegatePoint<2>>> {
        let results: Vec<DelegatePoint<2>> = unchecked_circles
            .par_iter()
            .filter_map(|circ| circ.get_possible_point(self.coord_tree(), self.bondlength()))
            .collect();
        if !results.is_empty() {
            Some(results)
        } else {
            None
        }
    }

    pub fn validate_site<T: Visualize + Clone>(&self, coord_site: T) -> Option<T> {
        let coord = coord_site.determine_coord();
        let bondlength = self.bondlength;
        let dist = bondlength.powi(2);
        if self
            .coord_tree
            .within_radius(&coord, bondlength)
            .iter()
            .any(|&&nb| {
                let distance = distance_squared(&coord, self.coord_tree.item(nb));
                matches!(approx_cmp_f64(distance, dist), FloatOrdering::Less)
            })
        {
            None
        } else {
            Some(coord_site)
        }
    }

    /// The first round search, abstract into sphere-sphere intersection.
    /// If the sphere does not have possible intersecting neighbours, then
    /// return early as `CoordResult::Sphere`. Otherwise, Use `CoordResult::Various`
    /// to unify the possible `CoordPoint` and `CoordCircle` (cut and intersect of two spheres)
    pub fn sphere_check(&self) -> SphereCheckResult {
        let mut results: Vec<CoordResult> = self
            .to_check
            .iter()
            .map(
                // Use `CoordResult::Various` to unify points and circles
                |&(atom_id, p)| -> CoordResult {
                    sphere_check_fn(atom_id, p, self.coord_tree(), self.bondlength())
                },
            )
            .collect();
        let unchecked_circles: Vec<Vec<CoordCircle>> = results
            .iter()
            .filter_map(|res| res.try_pull_circles_from_various().ok())
            .collect();
        let points: Vec<Vec<MultiCoordPoint>> = results
            .iter_mut()
            .filter_map(|res| res.try_pull_single_points_from_various().ok())
            .collect();
        SphereCheckResult::new(points.concat(), unchecked_circles.concat())
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
        coord_sites: Vec<T>,
        search_config: &SearchConfig,
    ) -> Vec<T> {
        coord_sites
            .into_iter()
            .filter_map(|coord_site| search_config.validate_site(coord_site))
            .collect()
    }
}

fn brute_force(
    origin: Point3<f64>,
    dist: f64,
    kdtree: &KdIndexTree<Point3<f64>>,
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
            .within_radius(&p, dist + 10_f64 * f64::EPSILON)
            .iter()
            .any(|&&nb| {
                let distance = distance_squared(&p, kdtree.item(nb));
                matches!(approx_cmp_f64(distance, dist.powi(2)), FloatOrdering::Less)
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
