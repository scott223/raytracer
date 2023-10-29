use std::f64::consts::PI;

use rand::{Rng, rngs::SmallRng};

use crate::{vec3::Vec3, onb::Onb, elements::{Hittable, Element}};

// struct for Probability Density Functions 
#[allow(dead_code)]
pub enum Pdf {
    SpherePDF(SpherePDF),
    CosinePDF(CosinePDF),
    HittablePDF(HittablePDF),
    MixedPDF(MixedPDF),
}

pub trait PDFTrait {
    fn generate(&self, rng: &mut SmallRng) -> Vec3;
    fn value(&self, direction: Vec3) -> f64;
}

impl PDFTrait for Pdf {
    fn generate(&self, rng: &mut SmallRng) -> Vec3 {
        match self {
            Pdf::SpherePDF(s) => s.generate(rng),
            Pdf::CosinePDF(c) => c.generate(rng),
            Pdf::HittablePDF(h) => h.generate(rng),
            Pdf::MixedPDF(m) => m.generate(rng),
        }
    }

    fn value(&self, direction: Vec3) -> f64 {
        match self {
            Pdf::SpherePDF(s) => s.value(direction),
            Pdf::CosinePDF(c) => c.value(direction),
            Pdf::HittablePDF(h) => h.value(direction),
            Pdf::MixedPDF(m) => m.value(direction),
        }
    }
}

pub struct MixedPDF {
    pub origin: Vec3,
    pub p1: Box<dyn PDFTrait>,
    pub p2: Box<dyn PDFTrait>,
}

// PDF for a hittable object
impl PDFTrait for MixedPDF {
    //combine the values from the two PDF's
    fn value(&self, direction: Vec3) -> f64 {
        0.5 * self.p1.value(direction) + 0.5 * self.p2.value(direction)
    }
    
    //pick a random ray from one of the PDFs
    fn generate(&self, rng: &mut SmallRng) -> Vec3 {
        let r = rng.gen_range(0.0..1.0);
        if r < 0.5 {
            self.p1.generate(rng)
        } else {
            self.p2.generate(rng)
        }
    }
}

impl MixedPDF {
    pub fn new(origin: Vec3, p1: impl PDFTrait + 'static, p2: impl PDFTrait + 'static) -> Self {
        MixedPDF { 
            origin,
            p1: Box::new(p1),
            p2: Box::new(p2),
         }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct HittablePDF {
    pub value: f64,
    pub origin: Vec3,
    pub object: Element,
}

// PDF for a hittable object
impl PDFTrait for HittablePDF {
    fn value(&self, direction: Vec3) -> f64 {
        self.object.pdf_value(self.origin, direction)
    }
    
    fn generate(&self, rng: &mut SmallRng) -> Vec3 {
        self.object.random(self.origin, rng)
    }
}

impl HittablePDF {
    pub fn new(origin: Vec3, object: Element) -> Self {
        HittablePDF { value: 0., 
            origin, 
            object,
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
    
    fn generate(&self, rng: &mut SmallRng) -> Vec3 {
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

    fn generate(&self, rng: &mut SmallRng) -> Vec3 {
        self.uvw.local_vec(Vec3::random_cosine_direction(rng))
    }
}