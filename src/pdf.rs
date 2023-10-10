use std::f64::consts::PI;

use rand::Rng;

use crate::{vec3::Vec3, onb::Onb};

// struct for Probability Density Functions 
#[derive(Debug, Clone, Copy)]
pub struct Pdf {
 //   pub value: Vec3,
}

pub trait PdfTrait {
    fn generate(&self, rng: &mut impl Rng) -> Vec3;
    fn value(&self, direction: Vec3) -> f64;
}

pub struct SpherePDF {
    pub value: f64,
}

// uniform density over a sphere
impl PdfTrait for SpherePDF {
    fn value(&self, _direction: Vec3) -> f64 {
        1. / (4. * PI)
    }
    
    fn generate(&self, rng: &mut impl Rng) -> Vec3 {
        Vec3::new_random_unit_vector(rng)
    }
}

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
impl PdfTrait for CosinePDF {
    fn value(&self, direction: Vec3) -> f64 {
        let cosine_theta = direction.normalized().dot(&self.uvw.w());
        0.0_f64.max(cosine_theta / PI)
    }

    fn generate(&self, rng: &mut impl Rng) -> Vec3 {
        self.uvw.local_vec(Vec3::random_cosine_direction(rng))
    }
}