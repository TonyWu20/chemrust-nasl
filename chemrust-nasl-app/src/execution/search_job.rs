use chemrust_geometry::Structure;
use chemrust_nasl::{SearchConfig, SearchReports};
use nalgebra::Point3;

use crate::error::RunError;
use crate::supportive_data::FractionalCoordRange;

use super::helpers::get_inbound_cartesian_coord;

#[derive(Debug, Clone, Copy, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub enum Axis {
    X,
    Y,
    Z,
}

pub trait SearchJob {
    fn target_bondlength(&self) -> f64;
    fn search_range(&self, axis: Axis) -> FractionalCoordRange;
    fn get_to_check_atom(&self, model: &Structure) -> Vec<(usize, Point3<f64>)> {
        let x_range = self.search_range(Axis::X);
        let y_range = self.search_range(Axis::Y);
        let z_range = self.search_range(Axis::Z);
        model
            .frac_coords
            .iter()
            .enumerate()
            .filter_map(|(i, coord)| {
                get_inbound_cartesian_coord(
                    coord,
                    x_range,
                    y_range,
                    z_range,
                    model.cell.as_ref().expect("the model carries a lattice"),
                )
                .map(|p| (i, p))
            })
            .collect()
    }
    fn all_atoms(&self, model: &Structure) -> Vec<Point3<f64>> {
        let all_range = FractionalCoordRange::new(0.0, 1.0);
        model
            .frac_coords
            .iter()
            .filter_map(|coord| {
                get_inbound_cartesian_coord(
                    coord,
                    all_range,
                    all_range,
                    all_range,
                    model.cell.as_ref().expect("the model carries a lattice"),
                )
            })
            .collect()
    }
    fn search(&self, model: &Structure) -> Result<SearchReports, RunError> {
        let to_check = self.get_to_check_atom(model);
        let all_atoms = self.all_atoms(model);
        let search_config = SearchConfig::new(&to_check, &all_atoms, self.target_bondlength());
        let search_report = search_config.search_sites();
        if search_report.viable_single_points().is_none()
            && search_report.viable_double_points().is_none()
            && search_report.points().is_none()
        {
            return Err(RunError::Message(
                "No available results for this config.".to_string(),
            ));
        }
        let validated_single_points = search_report
            .viable_single_points()
            .cloned()
            .map(|v| SearchReports::validated_results(v, &search_config));
        let validated_double_points = search_report
            .viable_double_points()
            .cloned()
            .map(|v| SearchReports::validated_results(v, &search_config));
        let validated_multi_points = search_report
            .points()
            .cloned()
            .map(|v| SearchReports::validated_results(v, &search_config));
        Ok(SearchReports::new(
            validated_multi_points,
            validated_single_points,
            validated_double_points,
        ))
    }
}
