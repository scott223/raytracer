use rand::rngs::SmallRng;
use serde::Deserialize;
use serde::Serialize;

use crate::bvh::Aabb;
use crate::linalg::Vec3;
use crate::materials::Material;
use crate::render::{Ray, Interval};

use super::Element;
use super::HitRecord;
use super::Hittable;

//Triangle element
#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct JSONTriangle {
    pub v0: Vec3,
    pub v1: Vec3,
    pub v2: Vec3,
    pub material: Material,
}

#[derive(Debug, Clone, Copy)]
pub struct Triangle {
    pub v0: Vec3,
    pub v1: Vec3,
    pub v2: Vec3,
    pub v0v1: Vec3,
    pub v0v2: Vec3,
    pub normal: Vec3,
    pub material: Material,
    pub bbox: Aabb,
}

impl JSONTriangle {
    pub fn add_as_element(&self, objects: &mut Vec<Element>) {
        let triangle = Element::Triangle(Triangle::new(self.v0, self.v1, self.v2, self.material));

        objects.push(triangle);
    }
}

impl Triangle {
    pub fn new_from_json_object(json_triangle: JSONTriangle) -> Self {
        Triangle::new(
            json_triangle.v0,
            json_triangle.v1,
            json_triangle.v2,
            json_triangle.material,
        )
    }

    // find the minimum xyz and max xyz coordinates and use for bounding box. add some padding as triangles are usually flat..
    // https://stackoverflow.com/questions/39974191/triangle-bounding-box
    pub fn new(v0: Vec3, v1: Vec3, v2: Vec3, material: Material) -> Self {
        let v0v1: Vec3 = v1 - v0;
        let v0v2: Vec3 = v2 - v0;
        let normal = v0v1.cross(&v0v2).normalized();

        Triangle {
            v0,
            v1,
            v2,
            v0v1,
            v0v2,
            normal,
            material,
            bbox: Aabb::new_from_points(v0.min(v1.min(v2)), v0.max(v1.max(v2))).pad(),
        }
    }
}

// Möller-Trumbore algorithm
// https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/moller-trumbore-ray-triangle-intersection.html
impl Hittable for Triangle {
    fn hit(&self, ray: &Ray, ray_t: &mut Interval) -> Option<HitRecord> {
        let pvec = ray.direction.cross(&self.v0v2);
        let det = self.v0v1.dot(&pvec);

        // applies culling
        if det < f64::EPSILON {
            return None;
        }

        let inv_det = 1. / det;

        let tvec = ray.origin - self.v0;
        let u = tvec.dot(&pvec) * inv_det;
        if !(0.0..=1.0).contains(&u) {
            return None;
        };

        let qvec = tvec.cross(&self.v0v1);
        let v = ray.direction.dot(&qvec) * inv_det;
        if v < 0.0 || u + v > 1.0 {
            return None;
        };

        let t = self.v0v2.dot(&qvec) * inv_det;

        // check if t is inside the interval (before the camera, and closer than an earlier hit)
        if !ray_t.contains(t) {
            return None;
        }

        // compute the intersection point
        let point: Vec3 = ray.at(t);

        // we have a hit, so we return a hit record
        Some(HitRecord {
            t,
            normal: self.normal,
            point,
            front_face: ray.direction.dot(&self.normal) < 0.0,
            material: self.material,
        })
    }

    fn bounding_box(&self) -> Aabb {
        self.bbox
    }

    fn pdf_value(&self, _origin: Vec3, _direction: Vec3) -> f64 {
        0.0
    }

    fn random(&self, _origin: Vec3, _rng: &mut SmallRng) -> Vec3 {
        Vec3::new(1.0, 0.0, 0.0)
    }
}
