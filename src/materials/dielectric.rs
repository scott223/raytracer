use rand::Rng;
use serde::{Deserialize, Serialize};

use crate::{color::Color, render::Ray, elements::HitRecord, linalg::Vec3};

use super::{Refracts, RefractRecord};

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct Dielectric {
    pub index_of_refraction: f64,
    // there is no albedo defined, as a glass material does not absorb color
    // albedo/atentuation is set to 1.0 1.0 1.0 in the code later
}

impl Dielectric {
    pub fn new(index_of_refraction: f64) -> Dielectric {
        Dielectric {
            index_of_refraction,
        }
    }
}

// source: Ray tracing in one Weekend
impl Refracts for Dielectric {
    fn refract(
        &self,
        ray: &Ray,
        hit_record: &HitRecord,
        rng: &mut impl Rng,
    ) -> Option<RefractRecord> {
        let albedo: Color = Color::new(1.0, 1.0, 1.0); // a glass material does not absorb any color/light so the albedo is 1.0
        let refraction_ratio: f64 = if hit_record.front_face {
            1.0 / self.index_of_refraction
        } else {
            self.index_of_refraction
        };
        let unit_direction: Vec3 = ray.direction.normalized(); // this should already be normalized, so we could remove this .normalize

        let minus_unit_direction: Vec3 = unit_direction * -1.0;
        let cos_theta = hit_record.normal.dot(&minus_unit_direction).min(1.0);
        let sin_theta = (1.0 - cos_theta * cos_theta).sqrt();

        let cannot_refract: bool = refraction_ratio * sin_theta > 1.0;

        if cannot_refract || super::metal::reflectance(cos_theta, refraction_ratio) > rng.gen::<f64>() {
            let direction: Vec3 = super::metal::reflect_vector(&unit_direction, &hit_record.normal);
            let reflected_ray: Ray = Ray::new(hit_record.point, direction);
            let refract_record: RefractRecord = RefractRecord {
                attenuation: albedo,
                ray: reflected_ray,
            };
            return Some(refract_record);
        } else {
            let direction: Vec3 =
                super::metal::refract_ray(&unit_direction, &hit_record.normal, refraction_ratio);
            let refracted_ray: Ray = Ray::new(hit_record.point, direction);
            let refract_record: RefractRecord = RefractRecord {
                attenuation: albedo,
                ray: refracted_ray,
            };
            return Some(refract_record);
        }
    }
}