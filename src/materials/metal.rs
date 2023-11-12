use rand::Rng;
use serde::{Deserialize, Serialize};

use crate::{color::Color, render::Ray, elements::HitRecord, linalg::Vec3};

use super::{Reflects, ReflectRecord};

// Metal material, with a fuzz factor. Metal reflects all rays in a predictable way (normal reflection)
#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct Metal {
    pub albedo: Color,
    pub fuzz: f64,
}

impl Metal {
    pub fn new(albedo: Color, fuzz: f64) -> Metal {
        Metal { albedo, fuzz }
    }
}

impl Reflects for Metal {
    fn reflect(
        &self,
        ray: &Ray,
        hit_record: &HitRecord,
        rng: &mut impl Rng,
    ) -> Option<ReflectRecord> {
        // get the direction of the reflected ray, and add a fuzz factor * a random unit vector
        // todo check if fuzz is between 0.0..1.0

        let new_direction = reflect_vector(&ray.direction.normalized(), &hit_record.normal)
            + Vec3::new_random_unit_vector(rng) * self.fuzz;

        if hit_record.normal.dot(&new_direction) > 0.0 {
            // the reflected ray, including fuzz unit sphere, is outside the material, so return a reflected ray
            let reflected = Ray::new(hit_record.point, new_direction);
            let reflect_record = ReflectRecord {
                attenuation: self.albedo,
                ray: Some(reflected),
            };

            return Some(reflect_record);
        } else {
            // return no ray, as the ray is absorbed by the material (due to fuzz factor)
            let reflect_record = ReflectRecord {
                attenuation: Color::new(0.0, 0.0, 0.0),
                ray: None,
            };

            return Some(reflect_record);
        }
    }
}

// reflect a ray with the same outgoing angle as incoming angle with the normal
pub fn reflect_vector(v: &Vec3, n: &Vec3) -> Vec3 {
    *v - *n * (2.0 * v.dot(n))
}

// refract a ray (for dielectric / glass materials)
// source: raytracing in one weekend
pub fn refract_ray(uv: &Vec3, n: &Vec3, etai_over_etat: f64) -> Vec3 {
    let minus_uv: Vec3 = *uv * -1.0;
    let cos_theta = n.dot(&minus_uv).min(1.0);
    let r_out_perp: Vec3 = (*uv + *n * cos_theta) * etai_over_etat;
    let r_out_parallel: Vec3 = *n * (((1.0 - r_out_perp.length_squared()).abs()).sqrt() * -1.0);
    r_out_perp + r_out_parallel
}

// check reflectance
// source: raytracing in one weekend
pub fn reflectance(cosine: f64, ref_idx: f64) -> f64 {
    let mut r0 = (1.0 - ref_idx) / (1.0 + ref_idx);
    r0 = r0 * r0;
    r0 + (1.0 - r0) * (1.0 - cosine).powi(5)
}