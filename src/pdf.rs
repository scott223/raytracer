use std::f64::consts::PI;

use rand::{rngs::SmallRng, Rng, seq::SliceRandom};

use crate::{
    elements::{Element, Hittable},
    onb::Onb,
    vec3::Vec3,
};

// struct for Probability Density Functions
#[allow(dead_code)]
pub enum Pdf<'a> {
    SpherePDF(SpherePDF),
    CosinePDF(CosinePDF),
    HittablePDF(HittablePDF<'a>),
    MixedPDF(MixedPDF<'a>),
}

pub trait PDFTrait {
    fn generate(&self, rng: &mut SmallRng) -> Vec3;
    fn value(&self, direction: Vec3) -> f64;
}

impl PDFTrait for Pdf<'_> {
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

// contains two PDFs
// need to specify lifetime so that references dont outlive the overal struct
pub struct MixedPDF<'a> {
    pub origin: Vec3,
    pub ratio: f64,
    pub p1: &'a Pdf<'a>,
    pub p2: &'a Pdf<'a>,
}

// PDF for a hittable object
impl PDFTrait for MixedPDF<'_> {
    //combine the values from the two PDF's
    fn value(&self, direction: Vec3) -> f64 {
        self.ratio * self.p1.value(direction) + (1.-self.ratio) * self.p2.value(direction)
    }

    //pick a random ray from one of the PDFs
    fn generate(&self, rng: &mut SmallRng) -> Vec3 {
        let r = rng.gen_range(0.0..1.0);
        if r < self.ratio {
            self.p1.generate(rng)
        } else {
            self.p2.generate(rng)
        }
    }
}

impl MixedPDF<'_> {
    // the lifetime of the mixed pdf p1 and p2 needs to be the same as the lifetime of the mixed pdf struct itself
    pub fn new<'a>(origin: Vec3, mut ratio: f64, p1: &'a Pdf, p2: &'a Pdf) -> MixedPDF<'a> {
        if !(0.0..=1.0).contains(&ratio) {
            ratio = 0.5;
        }

        MixedPDF {
            origin,
            ratio,
            p1: p1,
            p2: p2,
        }
        
    }
}

#[derive(Debug, Clone)]
pub struct HittablePDF<'a> {
    pub value: f64,
    pub origin: Vec3,
    pub objects: &'a Vec<&'a Element>,
}

// PDF for a hittable object
impl PDFTrait for HittablePDF<'_> {

    // returns the PDF value as weighted sum of all the values of the objects referenced
    fn value(&self, direction: Vec3) -> f64 {
        let weight: f64 = 1.0 / self.objects.len() as f64;
        let mut sum: f64 = 0.0;
        
        for o in self.objects.iter() {
            sum += weight * o.pdf_value(self.origin, direction);
        }

        sum
    }

    // generates a point on one of the (randomly picked) objects in the PDF
    fn generate(&self, rng: &mut SmallRng) -> Vec3 {
        self.objects.choose(rng).expect("No object returned with random point method").random(self.origin, rng)
    }
}

impl HittablePDF<'_> {
    pub fn new<'a>(origin: Vec3, objects: &'a Vec<&'a Element>) -> HittablePDF<'a> {
        HittablePDF {
            value: 0.,
            origin,
            objects,
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
