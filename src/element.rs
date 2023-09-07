use crate::materials::Material;
use crate::sphere::Sphere;
use crate::plane::Plane;
use crate::vec3::Vec3;
use crate::ray::Ray;

use serde::{Serialize, Deserialize};

// hitrecord gets returned on a hit, containting the point on the ray, the point in the global coordinate system, the normal and the material for the hit
#[derive(Debug)]
pub struct HitRecord {
    pub t: f64,
    pub point: Vec3,
    pub normal: Vec3,
    pub material: Material,
}

// hittable trait defines the hit function for each element
pub trait Hittable {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord>;
}

// enum for all the different elements
#[derive(Serialize, Deserialize, Debug)]
pub enum Element{
    Sphere(Sphere),
    Plane(Plane),
}

// matching the hit function with the hittable trait for each type of element
impl Hittable for Element {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord>{
        match *self {
            Element::Sphere(ref s) => s.hit(ray, t_min, t_max),
            Element::Plane(ref p) => p.hit(ray, t_min, t_max),
        }
    }
}