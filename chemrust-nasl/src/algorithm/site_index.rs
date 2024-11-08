use chemrust_core::data::atom::CoreAtomData;
use kiddo::ImmutableKdTree;
use nalgebra::Point3;

fn build_kd_tree_from_points(points: &[Point3<f64>]) -> ImmutableKdTree<f64, 3> {
    let entries = points.iter().map(|&p| p.into()).collect::<Vec<[f64; 3]>>();
    ImmutableKdTree::new_from_slice(&entries)
}

#[derive(Debug, Clone)]
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

impl<T: CoreAtomData> From<&T> for SiteIndex {
    fn from(value: &T) -> Self {
        let (coords, coords_entries): (Vec<Point3<f64>>, Vec<[f64; 3]>) = value
            .coords_repr()
            .iter()
            .map(|coord| -> (Point3<f64>, [f64; 3]) { (coord.raw_data(), coord.raw_data().into()) })
            .unzip();
        let coord_tree = ImmutableKdTree::new_from_slice(&coords_entries);
        SiteIndex { coords, coord_tree }
    }
}
