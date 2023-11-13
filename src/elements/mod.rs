mod triangle;
pub use triangle::JSONTriangle;
pub use triangle::Triangle;

mod sphere;
pub use sphere::JSONSphere;
pub use sphere::Sphere;

mod quad;
pub use quad::JSONQuad;
pub use quad::Quad;

mod obj;
pub use obj::JSONObj;

mod cube;
pub use cube::JSONCube;

use rand::rngs::SmallRng;
use serde::Deserialize;
use serde::Serialize;

use crate::bvh::Aabb;
use crate::linalg::Vec3;
use crate::materials::Material;
use crate::render::Interval;
use crate::render::Ray;

// hitrecord gets returned on a hit, containting the point on the ray, the point in the global coordinate system, the normal and the material for the hit
#[derive(Debug)]
pub struct HitRecord {
    pub t: f64,
    pub point: Vec3,
    pub normal: Vec3,
    pub front_face: bool,
    pub material: Material,
}

// hittable trait defines the hit function for each element
pub trait Hittable {
    fn hit(&self, ray: &Ray, ray_t: &mut Interval) -> Option<HitRecord>;
    fn bounding_box(&self) -> Aabb;
    fn pdf_value(&self, origin: Vec3, direction: Vec3) -> f64;
    fn random(&self, origin: Vec3, rng: &mut SmallRng) -> Vec3;
}

// enum for all the different elements, simplified JSON representation
#[derive(Serialize, Deserialize, Debug, Clone)]
pub enum JSONElement {
    JSONSphere(JSONSphere),
    JSONQuad(JSONQuad),
    JSONTriangle(JSONTriangle),
    JSONBox(JSONCube),
    JSONObj(JSONObj),
}

// enum for all the different elements
#[derive(Debug, Clone, Copy)]
pub enum Element {
    Sphere(Sphere),
    Quad(Quad),
    Triangle(Triangle),
}

//implement the methods for the Element
impl Element {
    // returns true if marked as an attractor
    pub fn is_attractor(&self) -> bool {
        match *self {
            Element::Sphere(ref s) => s.is_attractor(),
            Element::Quad(ref q) => q.is_attractor(),
            _ => false,
        }
    }
}

// element transformations
#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct Transpose {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct Rotate {
    pub theta_x: f64,
    pub theta_y: f64,
    pub theta_z: f64,
}

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct Scale {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

// matching the hit function with the hittable trait for each type of element
impl Hittable for Element {
    fn hit(&self, ray: &Ray, ray_t: &mut Interval) -> Option<HitRecord> {
        match *self {
            Element::Sphere(ref s) => s.hit(ray, ray_t),
            Element::Quad(ref q) => q.hit(ray, ray_t),
            Element::Triangle(ref t) => t.hit(ray, ray_t),
        }
    }

    fn bounding_box(&self) -> Aabb {
        match *self {
            Element::Sphere(ref s) => s.bounding_box(),
            Element::Quad(ref q) => q.bounding_box(),
            Element::Triangle(ref t) => t.bounding_box(),
        }
    }

    fn pdf_value(&self, origin: Vec3, direction: Vec3) -> f64 {
        match *self {
            Element::Sphere(ref s) => s.pdf_value(origin, direction),
            Element::Quad(ref q) => q.pdf_value(origin, direction),
            Element::Triangle(ref t) => t.pdf_value(origin, direction),
        }
    }

    fn random(&self, origin: Vec3, rng: &mut SmallRng) -> Vec3 {
        match *self {
            Element::Sphere(ref s) => s.random(origin, rng),
            Element::Quad(ref q) => q.random(origin, rng),
            Element::Triangle(ref t) => t.random(origin, rng),
        }
    }
}
