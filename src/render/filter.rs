/// Enum for different filter classes
#[derive(Debug, Clone)]
pub enum Filter {
    MitchellNetravali(MitchellNetravali),
}

impl Filter {
    pub fn evaluate(&self, sample_x: f64, sample_y: f64) -> f64 {
        match self {
            Filter::MitchellNetravali(filter) => filter.evaluate(sample_x, sample_y),
        }
    }
    pub fn get_radius(&self) -> (f64, f64) {
        match self {
            Filter::MitchellNetravali(filter) => filter.get_radius(),
        }
    }
}

/// Struct for MitchellNetravali filter
#[derive(Debug, Default, Copy, Clone)]
pub struct MitchellNetravali {
    width: f64,
    height: f64,
    inv_width: f64,
    inv_height: f64,
    b: f64,
    c: f64,
}

impl MitchellNetravali {
    /// Returns a new MN filter
    pub fn new(w: f64, h: f64, b: f64, c: f64) -> Self {
        MitchellNetravali {
            width: w,
            height: h,
            inv_width: 1.0 / w,
            inv_height: 1.0 / h,
            b,
            c,
        }
    }

    #[inline(always)]
    fn mitchell_1d(&self, x: f64) -> f64 {
        let fx = x.abs() * 2.0;
        if fx < 1.0 {
            ((12.0 - 9.0 * self.b - 6.0 * self.c) * fx * fx * fx
                + (-18.0 + 12.0 * self.b + 6.0 * self.c) * fx * fx
                + (6.0 - 2.0 * self.b))
                * (1.0 / 6.0)
        } else if fx < 2.0 {
            ((-self.b - 6.0 * self.c) * fx * fx * fx
                + (6.0 * self.b + 30.0 * self.c) * fx * fx
                + (-12.0 * self.b - 48.0 * self.c) * fx
                + (8.0 * self.b + 24.0 * self.c))
                * (1.0 / 6.0)
        } else {
            0.0
        }
    }

    /// Evaluates the MN filter for a given sample position
    #[inline(always)]
    pub fn evaluate(&self, sample_x: f64, sample_y: f64) -> f64 {
        self.mitchell_1d(sample_x * self.inv_width) * self.mitchell_1d(sample_y * self.inv_height)
    }

    /// Returns the radius set for the MN filter
    #[inline(always)]
    pub fn get_radius(&self) -> (f64, f64) {
        (self.width, self.height)
    }
}
