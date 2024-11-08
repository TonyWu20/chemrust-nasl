use nalgebra::Point3;

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
