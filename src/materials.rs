use crate::color::Color;
use crate::elements::HitRecord;
use crate::ray::Ray;
use crate::vec3::Vec3;
use rand::Rng;

use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub enum Material {
    Lambertian(Lambertian),
    Metal(Metal),
    Dielectric(Dielectric),
    DiffuseLight(DiffuseLight),
}

// trait for a material that scatters
pub trait Scatterable {
    fn scatter(&self, ray: &Ray, hit_record: &HitRecord) -> Option<(Option<Ray>, Color)>;
}

pub trait Emmits {
    fn emitted(&self, ray: &Ray, hit_record: &HitRecord) -> Option<Color>;
}

// link the trait implementation to the materials
// we now assume every material scatters, so each material needs a scatter function
impl Scatterable for Material {
    fn scatter(&self, ray: &Ray, hit_record: &HitRecord) -> Option<(Option<Ray>, Color)> {
        match self {
            Material::Lambertian(l) => l.scatter(ray, hit_record),
            Material::Metal(m) => m.scatter(ray, hit_record),
            Material::Dielectric(d) => d.scatter(ray, hit_record),
            _ => None,
        }
    }
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

// create a new Lambertian material
impl DiffuseLight {
    pub fn new(albedo: Color) -> DiffuseLight {
        DiffuseLight { albedo }
    }
}

impl Emmits for DiffuseLight {
    fn emitted(&self, _ray: &Ray, _hit_record: &HitRecord) -> Option<Color> {
        Some(self.albedo)
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
    fn scatter(&self, _ray: &Ray, hit_record: &HitRecord) -> Option<(Option<Ray>, Color)> {
        let mut new_direction = hit_record.normal + Vec3::new_random_unit_vector(); //lambertian distribution

        // if the direction is almost zero, scatter to the normal
        if new_direction.near_zero() {
            new_direction = hit_record.normal;
        }

        // create the new ray
        let scattered = Ray::new(hit_record.point, new_direction);

        Some((Some(scattered), self.albedo))
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

// reflect a ray with the same outgoing angle as incoming angle with the normal
fn reflect(v: &Vec3, n: &Vec3) -> Vec3 {
    *v - *n * (2.0 * v.dot(n))
}

impl Scatterable for Metal {
    // create a reflected ray
    fn scatter(&self, ray: &Ray, hit_record: &HitRecord) -> Option<(Option<Ray>, Color)> {
        // get the direction of the reflected ray, and add a fuzz factor * a random unit vector
        let new_direction = reflect(&ray.direction.normalized(), &hit_record.normal)
            + Vec3::new_random_unit_vector() * self.fuzz;

        if hit_record.normal.dot(&new_direction) > 0.0 {
            // the reflected ray, including fuzz unit sphere, is outside the material, so return a reflected ray
            let reflected = Ray::new(hit_record.point, new_direction);
            Some((Some(reflected), self.albedo))
        } else {
            // return no ray, as the ray is absorbed by the material (due to fuzz factor)
            Some((None, self.albedo))
        }
    }
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
    fn scatter(&self, ray: &Ray, hit_record: &HitRecord) -> Option<(Option<Ray>, Color)> {
        let mut rng = rand::thread_rng();
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
            Some((Some(reflected_ray), albedo))
        } else {
            let direction: Vec3 = refract(&unit_direction, &hit_record.normal, refraction_ratio);
            let refracted_ray: Ray = Ray::new(hit_record.point, direction);
            Some((Some(refracted_ray), albedo))
        }
    }
}


