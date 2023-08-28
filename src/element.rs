use crate::materials::Material;
use crate::sphere::Sphere;
use crate::plane::Plane;
use crate::vec3::Vec3;
use crate::ray::Ray;

#[derive(Debug)]
pub struct HitRecord {
    pub t: f64,
    pub point: Vec3,
    pub normal: Vec3,
    pub material: Material,
}

pub trait Hittable {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord>;
}

#[derive(Debug)]
pub enum Element{
    Sphere(Sphere),
    Plane(Plane),
}

impl Hittable for Element {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord>{
        match *self {
            Element::Sphere(ref s) => s.hit(ray, t_min, t_max),
            Element::Plane(ref p) => p.hit(ray, t_min, t_max),
        }
    }

}