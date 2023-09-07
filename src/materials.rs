use crate::color::Color;
use crate::element::HitRecord;
use crate::vec3::Vec3;
use crate::ray::Ray;

use serde::{Serialize, Deserialize};

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub enum Material {
    Lambertian(Lambertian),
    Metal(Metal),
    Dielectric(Dielectric),
}

// trait for a material that scatters
pub trait Scatterable {
    fn scatter(&self, ray: &Ray, hit_record: &HitRecord) -> Option<(Option<Ray>, Color)>;
}

// linkt the trait implementation to the materials
impl Scatterable for Material {
    fn scatter(&self, ray: &Ray, hit_record: &HitRecord) -> Option<(Option<Ray>, Color)> {
        match self {
            Material::Lambertian(l) => l.scatter(ray, hit_record),
            Material::Metal(m) => m.scatter(ray, hit_record),
            Material::Dielectric(d) => d.scatter(ray, hit_record), 
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
    fn scatter(&self, _ray: &Ray, hit_record: &HitRecord) -> Option<(Option<Ray>, Color)> {
        let mut new_direction = hit_record.normal + Vec3::new_random_unit_vector(); //lambertian distribution
        // if the direction is almost zero, scatter to the normal
        if new_direction.near_zero() {
            new_direction = hit_record.normal;
        }

        let scattered = Ray::new(hit_record.point, new_direction);
        let albedo = self.albedo;

        Some((Some(scattered), albedo))
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
        Metal { 
            albedo,
            fuzz, 
        }
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
        let new_direction = reflect(&ray.direction.normalized(), &hit_record.normal) + Vec3::new_random_unit_vector() * self.fuzz;
        let albedo = self.albedo;

        if hit_record.normal.dot(&new_direction) > 0.0 { 
            // the reflected ray, including fuzz unit sphere, is outside the material, so return a reflected ray
            let reflected = Ray::new(hit_record.point, new_direction);
            Some((Some(reflected), albedo))
        } else {
            // return no ray, as the ray is absorbed by the material (due to fuzz factor)
            Some((None, albedo)) 
        }
    }
}

    // refract a ray (for dielectric / glass materials)
fn refract(uv: &Vec3, n: &Vec3, etai_over_etat: f64) -> Vec3 {
    let minus_uv: Vec3 = *uv * -1.0;
    let cos_theta = n.dot(&minus_uv).min(1.0);
    let r_out_perp: Vec3 = (*uv + *n * cos_theta) * etai_over_etat;
    let r_out_parallel: Vec3 = *n * (((1.0 - r_out_perp.length_squared()).abs()).sqrt() * -1.0);
    r_out_perp + r_out_parallel
}

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct Dielectric {
    pub index_of_refraction: f64,
}

impl Dielectric {
    pub fn new(index_of_refraction: f64) -> Dielectric {
        Dielectric { 
            index_of_refraction,
        }
    }
}

impl Scatterable for Dielectric {

    fn scatter(&self, ray: &Ray, hit_record: &HitRecord) -> Option<(Option<Ray>, Color)> {
        let albedo: Color = Color::new(1.0, 1.0, 1.0); 
        let refraction_ratio: f64 = 1.0 / self.index_of_refraction;
        let unit_direction: Vec3 = ray.direction.normalized(); // this should already be normalized, so we could remove this .normalize

        let minus_unit_direction: Vec3 = unit_direction * -1.0;
        let cos_theta = hit_record.normal.dot(&minus_unit_direction).min(1.0);
        let sin_theta = (1.0 - cos_theta*cos_theta).sqrt();

        let cannot_refract: bool = refraction_ratio * sin_theta > 1.0;

        if cannot_refract {
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