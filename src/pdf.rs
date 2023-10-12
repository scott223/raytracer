use std::f64::consts::PI;

use rand::Rng;

use crate::{vec3::Vec3, onb::Onb};

// struct for Probability Density Functions 
#[derive(Debug, Clone, Copy)]
pub enum PDF {
    SpherePDF(SpherePDF),
    CosinePDF(CosinePDF),
}

pub trait PDFTrait {
    fn generate(&self, rng: &mut impl Rng) -> Vec3;
    fn value(&self, direction: Vec3) -> f64;
}

impl PDFTrait for PDF {
    fn generate(&self, rng: &mut impl Rng) -> Vec3 {
        match self {
            PDF::SpherePDF(s) => s.generate(rng),
            PDF::CosinePDF(c) => c.generate(rng),
        }
    }

    fn value(&self, direction: Vec3) -> f64 {
        match self {
            PDF::SpherePDF(s) => s.value(direction),
            PDF::CosinePDF(c) => c.value(direction),
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct SpherePDF {
    pub value: f64,
}

// uniform density over a sphere
impl PDFTrait for SpherePDF {
    fn value(&self, _direction: Vec3) -> f64 {
        1. / (4. * PI)
    }
    
    fn generate(&self, rng: &mut impl Rng) -> Vec3 {
        Vec3::new_random_unit_vector(rng)
    }
}

#[derive(Debug, Clone, Copy)]
pub struct CosinePDF {
    pub uvw: Onb,
}

impl CosinePDF {
    pub fn new(w: Vec3) -> Self {
        Self {
            uvw: Onb::build_from_w(w),
        }
    }
}

// cosine density
impl PDFTrait for CosinePDF {
    fn value(&self, direction: Vec3) -> f64 {
        let cosine_theta = direction.normalized().dot(&self.uvw.w());
        0.0_f64.max(cosine_theta / PI)
    }

    fn generate(&self, rng: &mut impl Rng) -> Vec3 {
        self.uvw.local_vec(Vec3::random_cosine_direction(rng))
    }
}