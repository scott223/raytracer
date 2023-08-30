use crate::color::Color;
use crate::element::HitRecord;
use crate::vec3::Vec3;
use crate::ray::Ray;

#[derive(Debug, Clone, Copy)]
pub enum Material {
    Lambertian(Lambertian),
    Metal(Metal),
}

pub trait Scatterable {
    fn scatter(&self, ray: &Ray, hit_record: &HitRecord) -> Option<(Option<Ray>, Color)>;
}

impl Scatterable for Material {
    fn scatter(&self, ray: &Ray, hit_record: &HitRecord) -> Option<(Option<Ray>, Color)> {
        match self {
            Material::Lambertian(l) => l.scatter(ray, hit_record),
            Material::Metal(m) => m.scatter(ray, hit_record),
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Lambertian {
    pub albedo: Color,
}

impl Lambertian {
    pub fn new(albedo: Color) -> Lambertian {
        Lambertian { albedo }
    }
}

impl Scatterable for Lambertian {
    fn scatter(&self, _ray: &Ray, hit_record: &HitRecord) -> Option<(Option<Ray>, Color)> {
        let mut new_direction = hit_record.normal + Vec3::new_random_unit_vector(); //lambertian distribution
        if new_direction.near_zero() {
            new_direction = hit_record.normal;
        }

        let scattered = Ray::new(hit_record.point, new_direction);
        let albedo = self.albedo;

        Some((Some(scattered), albedo))
    }
}

#[derive(Debug, Clone, Copy)]
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

fn reflect(v: &Vec3, n: &Vec3) -> Vec3 {
    *v - *n * (2.0 * v.dot(n))
}

impl Scatterable for Metal {
    fn scatter(&self, ray: &Ray, hit_record: &HitRecord) -> Option<(Option<Ray>, Color)> {
        // get the direction of the reflected ray, and add a fuzz factor * a random unit vector
        let new_direction = reflect(&ray.direction.normalized(), &hit_record.normal) + Vec3::new_random_unit_vector() * self.fuzz;
        let albedo = self.albedo;

        if hit_record.normal.dot(&new_direction) > 0.0 { 
            // the reflected ray, including fuzz unit sphere, is outside the material, so return a reflected ray
            let scattered = Ray::new(hit_record.point, new_direction);
            Some((Some(scattered), albedo))
        } else {
            // return no ray, as the ray is absorbed by the material (due to fuzz factor)
            Some((None, albedo)) 
        }
    }
}