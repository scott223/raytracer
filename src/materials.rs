use std::cell::Ref;
use std::f64::consts::PI;

use crate::elements::HitRecord;
use crate::pdf::{CosinePDF, Pdf};
use crate::ray::Ray;
use crate::vec3::Vec3;
use crate::{color::Color, onb::Onb};

use rand::Rng;
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub enum Material {
    Lambertian(Lambertian),
    Metal(Metal),
    Dielectric(Dielectric),
    DiffuseLight(DiffuseLight),
}

//need to specify lifetime for the pdf, as there might be a mixed pdf in there with a reference to other pdfs
pub struct ScatterRecord<'a> {
    pub attenuation: Color,
    pub pdf: Pdf<'a>,
}

// trait for a material that scatters
pub trait Scatterable {
    fn scatter(
        &self,
        ray: &Ray,
        hit_record: &HitRecord,
        rng: &mut impl Rng,
    ) -> Option<ScatterRecord>;
    //fn scattering_pdf(&self, ray: &Ray, hit_record: &HitRecord, scattered_ray: &Ray) -> f64;
}

// link the trait implementation to the materials
impl Scatterable for Material {
    fn scatter(
        &self,
        ray: &Ray,
        hit_record: &HitRecord,
        rng: &mut impl Rng,
    ) -> Option<ScatterRecord> {
        match self {
            Material::Lambertian(l) => l.scatter(ray, hit_record, rng),
            _ => None,
        }
    }
}

pub struct ReflectRecord {
    pub attenuation: Color,
    pub ray: Option<Ray>,
}

// trait for a material that reflects
pub trait Reflects {
    fn reflect(
        &self,
        ray: &Ray,
        hit_record: &HitRecord,
        rng: &mut impl Rng,
    ) -> Option<ReflectRecord>;
    //fn scattering_pdf(&self, ray: &Ray, hit_record: &HitRecord, scattered_ray: &Ray) -> f64;
}

// link the trait implementation to the materials
// we now assume every material scatters, so each material needs a scatter function
impl Reflects for Material {
    fn reflect(
        &self,
        ray: &Ray,
        hit_record: &HitRecord,
        rng: &mut impl Rng,
    ) -> Option<ReflectRecord> {
        match self {
            Material::Metal(m) => m.reflect(ray, hit_record, rng),
            _ => None,
        }
    }
}

pub trait Emmits {
    fn emitted(&self, ray: &Ray, hit_record: &HitRecord) -> Option<Color>;
}

impl Emmits for Material {
    fn emitted(&self, ray: &Ray, hit_record: &HitRecord) -> Option<Color> {
        match self {
            Material::DiffuseLight(d) => d.emitted(ray, hit_record),
            _ => None,
        }
    }
}

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct DiffuseLight {
    // albedo is defined as amount of color not absorbed
    pub albedo: Color,
}

// create a new Diffuse light material
impl DiffuseLight {
    pub fn new(albedo: Color) -> DiffuseLight {
        DiffuseLight { albedo }
    }
}

impl Emmits for DiffuseLight {
    fn emitted(&self, _ray: &Ray, hit_record: &HitRecord) -> Option<Color> {
        if hit_record.front_face {
            Some(self.albedo)
        } else {
            None
        }
    }
}

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
        rng: &mut impl Rng,
    ) -> Option<ScatterRecord> {
        //lambertian pdf

        let scatter: ScatterRecord = ScatterRecord {
            attenuation: self.albedo,
            pdf: Pdf::CosinePDF(CosinePDF::new(hit_record.normal)),
        };

        Some(scatter)
    }
}

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
                attenuation: self.albedo,
                ray: None,
            };

            return Some(reflect_record);
        }
    }
}

// reflect a ray with the same outgoing angle as incoming angle with the normal
fn reflect_vector(v: &Vec3, n: &Vec3) -> Vec3 {
    *v - *n * (2.0 * v.dot(n))
}

// refract a ray (for dielectric / glass materials)
// source: raytracing in one weekend
fn refract(uv: &Vec3, n: &Vec3, etai_over_etat: f64) -> Vec3 {
    let minus_uv: Vec3 = *uv * -1.0;
    let cos_theta = n.dot(&minus_uv).min(1.0);
    let r_out_perp: Vec3 = (*uv + *n * cos_theta) * etai_over_etat;
    let r_out_parallel: Vec3 = *n * (((1.0 - r_out_perp.length_squared()).abs()).sqrt() * -1.0);
    r_out_perp + r_out_parallel
}

// check reflectance
// source: raytracing in one weekend
fn reflectance(cosine: f64, ref_idx: f64) -> f64 {
    let mut r0 = (1.0 - ref_idx) / (1.0 + ref_idx);
    r0 = r0 * r0;
    r0 + (1.0 - r0) * (1.0 - cosine).powi(5)
}

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct Dielectric {
    pub index_of_refraction: f64,
    // there is no albedo defined, as a glass material does not absorb color
}

impl Dielectric {
    pub fn new(index_of_refraction: f64) -> Dielectric {
        Dielectric {
            index_of_refraction,
        }
    }
}

// source: Ray tracing in one Weekend
impl Scatterable for Dielectric {
    fn scatter(
        &self,
        ray: &Ray,
        hit_record: &HitRecord,
        rng: &mut impl Rng,
    ) -> Option<ScatterRecord> {
        //let mut rng = rand::thread_rng();
        /*
        let albedo: Color = Color::new(1.0, 1.0, 1.0); // a glass material does not absorb any color/light so the albedo is 1.0
        let refraction_ratio: f64 = if hit_record.front_face { 1.0 / self.index_of_refraction } else { self.index_of_refraction };
        let unit_direction: Vec3 = ray.direction.normalized(); // this should already be normalized, so we could remove this .normalize

        let minus_unit_direction: Vec3 = unit_direction * -1.0;
        let cos_theta = hit_record.normal.dot(&minus_unit_direction).min(1.0);
        let sin_theta = (1.0 - cos_theta * cos_theta).sqrt();

        let cannot_refract: bool = refraction_ratio * sin_theta > 1.0;

        if cannot_refract || reflectance(cos_theta, refraction_ratio) > rng.gen::<f64>() {
            let direction: Vec3 = reflect(&unit_direction, &hit_record.normal);
            let reflected_ray: Ray = Ray::new(hit_record.point, direction);
            Some((Some(reflected_ray), Some(1.), albedo))
        } else {
            let direction: Vec3 = refract(&unit_direction, &hit_record.normal, refraction_ratio);
            let refracted_ray: Ray = Ray::new(hit_record.point, direction);
            Some((Some(refracted_ray), Some(1.0), albedo))
        }
        */
        None
    }
}

//TODO Tests
