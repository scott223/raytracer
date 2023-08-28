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
    fn scatter(&self, ray: &Ray, hit_record: &HitRecord) -> Option<(Option<Ray>, Color)> {
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
}

impl Metal {
    pub fn new(albedo: Color) -> Metal {
        Metal { albedo }
    }
}

fn reflect(v: &Vec3, n: &Vec3) -> Vec3 {
    *v - *n * (2.0 * v.dot(n))
}

impl Scatterable for Metal {
    fn scatter(&self, ray: &Ray, hit_record: &HitRecord) -> Option<(Option<Ray>, Color)> {
        let new_direction = reflect(&ray.direction.normalized(), &hit_record.normal);

        let scattered = Ray::new(hit_record.point, new_direction);
        let albedo = self.albedo;

        Some((Some(scattered), albedo))
    }
}