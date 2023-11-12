use rand::Rng;
use serde::{Deserialize, Serialize};

use crate::{color::Color, render::{Ray, Pdf, CosinePDF}, elements::HitRecord};

use super::{Scatterable, ScatterRecord};

// Lambertian (diffuse) material, that scatters rays in a semi-random direction (lambertian distribution = more concentrated around the normal)
#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct Lambertian {
    // albedo is defined as amount of color not absorbed
    pub albedo: Color,
}

// create a new Lambertian material
impl Lambertian {
    pub fn new(albedo: Color) -> Lambertian {
        Lambertian { albedo }
    }
}

impl Scatterable for Lambertian {
    // create a scattered ray, randomized but with a lambartian distribution around the normal
    fn scatter(
        &self,
        _ray: &Ray,
        hit_record: &HitRecord,
        _rng: &mut impl Rng,
    ) -> Option<ScatterRecord> {
        //lambertian pdf

        let scatter: ScatterRecord = ScatterRecord {
            attenuation: self.albedo,
            pdf: Pdf::CosinePDF(CosinePDF::new(hit_record.normal)),
        };

        Some(scatter)
    }
}