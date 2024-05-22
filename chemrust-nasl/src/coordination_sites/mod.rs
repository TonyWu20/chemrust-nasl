mod coord_circle;
mod coord_point;
mod coord_sphere;
mod visualize;

pub use coord_circle::CoordCircle;
pub use coord_point::CoordPoint;
pub use coord_sphere::CoordSphere;
pub use visualize::*;

pub trait CoordSite {
    fn connecting_atoms_msg(&self) -> String;
}

impl CoordSite for CoordCircle {
    fn connecting_atoms_msg(&self) -> String {
        format!(
            "double_{}",
            self.atom_ids()
                .iter()
                .map(|v| format!("{v}"))
                .collect::<Vec<String>>()
                .join("_")
        )
    }
}

impl CoordSite for CoordSphere {
    fn connecting_atoms_msg(&self) -> String {
        format!("single_{}", self.atom_id)
    }
}

impl CoordSite for CoordPoint {
    fn connecting_atoms_msg(&self) -> String {
        format!(
            "multi_cn_{}_{}",
            self.atom_ids().len(),
            self.atom_ids()
                .iter()
                .map(|v| format!("{v}"))
                .collect::<Vec<String>>()
                .join("_")
        )
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum CoordResult {
    Invalid,
    Empty,
    Sphere(CoordSphere),
    Circle(CoordCircle),
    SinglePoint(CoordPoint),
    Points(Vec<CoordPoint>),
    Various(Vec<CoordResult>),
}

impl CoordResult {
    pub fn try_into_sphere(&self) -> Result<CoordSphere, &Self> {
        if let Self::Sphere(v) = self {
            Ok(*v)
        } else {
            Err(self)
        }
    }

    pub fn try_into_circle(self) -> Result<CoordCircle, Self> {
        if let Self::Circle(v) = self {
            Ok(v)
        } else {
            Err(self)
        }
    }

    pub fn try_into_single_point(self) -> Result<CoordPoint, Self> {
        if let Self::SinglePoint(v) = self {
            Ok(v)
        } else {
            Err(self)
        }
    }

    pub fn try_into_points(self) -> Result<Vec<CoordPoint>, Self> {
        if let Self::Points(v) = self {
            Ok(v)
        } else {
            Err(self)
        }
    }
    pub fn try_pull_single_points_from_various(&self) -> Result<Vec<CoordPoint>, &Self> {
        if let Self::Various(v) = self {
            Ok(v.iter()
                .filter_map(|res| {
                    if let CoordResult::SinglePoint(p) = res {
                        Some(p.to_owned())
                    } else {
                        None
                    }
                })
                .collect())
        } else {
            Err(self)
        }
    }
    pub fn try_pull_points_from_various(&self) -> Result<Vec<CoordPoint>, &Self> {
        if let Self::Various(v) = self {
            let points: Vec<Vec<CoordPoint>> = v
                .iter()
                .filter_map(|res| {
                    if let CoordResult::Points(p) = res {
                        Some(p.to_owned())
                    } else {
                        None
                    }
                })
                .collect();
            Ok(points.concat())
        } else {
            Err(self)
        }
    }
    pub fn try_pull_circles_from_various(&self) -> Result<Vec<CoordCircle>, &Self> {
        if let Self::Various(v) = self {
            let circles: Vec<CoordCircle> = v
                .iter()
                .filter_map(|res| {
                    if let CoordResult::Circle(c) = res {
                        Some(*c)
                    } else {
                        None
                    }
                })
                .collect();
            Ok(circles)
        } else {
            Err(self)
        }
    }
}
