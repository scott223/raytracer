use rand::rngs::SmallRng;
use rand::Rng;
use serde::Deserialize;
use serde::Serialize;

use crate::bvh::Aabb;
use crate::linalg::Vec3;
use crate::materials::Material;
use crate::render::{Interval, Ray};

use super::Element;
use super::HitRecord;
use super::Hittable;

//Quad element
// simplified JSON version of the quad
#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct JSONQuad {
    pub q: Vec3,
    pub u: Vec3,
    pub v: Vec3,
    pub material: Material,
    pub attractor: Option<bool>,
}

impl JSONQuad {
    pub fn add_as_element(&self, objects: &mut Vec<Element>) {
        let object = Element::Quad(Quad::new(
            self.q,
            self.u,
            self.v,
            self.material,
            self.attractor,
        ));
        objects.push(object);
    }
}

// actual quad element
#[derive(Debug, Clone, Copy)]
pub struct Quad {
    pub q: Vec3,
    pub u: Vec3,
    pub v: Vec3,
    pub material: Material,
    pub attractor: Option<bool>,

    //following fields are used for every hit calculation, so we pre-calculate in the constructor
    pub n: Vec3,
    pub normal: Vec3,
    pub d: f64,
    pub w: Vec3,
    pub area: f64,
    pub bbox: Aabb,
}

impl Quad {
    pub fn new_from_json_object(json_quad: JSONQuad) -> Self {
        Quad::new(
            json_quad.q,
            json_quad.u,
            json_quad.v,
            json_quad.material,
            json_quad.attractor,
        )
    }

    // Creating a new Quad with lower left point Q and vectors u and v
    pub fn new(q: Vec3, u: Vec3, v: Vec3, material: Material, attractor: Option<bool>) -> Self {
        let n: Vec3 = u.cross(&v);
        let normal: Vec3 = n.normalized();
        let d: f64 = normal.dot(&q);
        let w: Vec3 = n / n.dot(&n);
        let area: f64 = n.length();

        Quad {
            q,
            u,
            v,
            n,
            normal,
            d,
            w,
            area,
            material,
            bbox: Aabb::new_from_points(q, q + u + v).pad(),
            attractor,
        }
    }

    //given the hit point in plane coordinates, return none if it is outside the primitive
    pub fn is_interior(a: f64, b: f64) -> Option<Vec<f64>> {
        if !(0.0..=1.0).contains(&a) || !(0.0..=1.0).contains(&b) {
            return None;
        }

        // we have this work around as we dont determine the u and v in the hitrecord yet
        Some(vec![a, b])
    }

    // returns true if self.attractor is set to true
    pub fn is_attractor(&self) -> bool {
        matches!(self.attractor, Some(true))
    }
}

impl Hittable for Quad {
    fn hit(&self, ray: &Ray, ray_t: &mut Interval) -> Option<HitRecord> {
        let denom = self.normal.dot(&ray.direction);

        // no hits, as the ray is parallel to the plane
        if denom.abs() < f64::EPSILON {
            return None;
        }

        // return false if the hit point paramater t is outside the ray interval
        let t: f64 = (self.d - self.normal.dot(&ray.origin)) / denom;
        if !ray_t.contains(t) {
            return None;
        }

        let intersection = ray.at(t);
        let planar_hitpoint_vector: Vec3 = intersection - self.q;
        let alpha: f64 = self.w.dot(&planar_hitpoint_vector.cross(&self.v));
        let beta: f64 = self.w.dot(&self.u.cross(&planar_hitpoint_vector));

        match Quad::is_interior(alpha, beta) {
            Some(_result) => {
                // nothing needed here i think
            }
            _ => {
                //is not in the interior, so break
                return None;
            }
        }

        // we have a hit, so we return a hit record
        Some(HitRecord {
            t,
            normal: self.normal,
            point: intersection,
            front_face: ray.direction.dot(&self.normal) < 0.0,
            material: self.material,
        })
    }

    fn bounding_box(&self) -> Aabb {
        //Aabb::new_from_points(self.Q, self.Q + self.u + self.v).pad()

        //we now calculate the bbox on initialization, so we can just return the field
        self.bbox
    }

    // returns the value for the probability distribution function for a given origing and direction
    fn pdf_value(&self, origin: Vec3, direction: Vec3) -> f64 {
        // check if this ray actually hits this quad
        let r = Ray::new(origin, direction);
        if let Some(hit) = self.hit(&r, &mut Interval::new(0.001, f64::INFINITY)) {
            let distance_squared = hit.t * hit.t * direction.length_squared();
            let cosine = direction.dot(&hit.normal).abs() / direction.length();

            distance_squared / (cosine * self.area)
        } else {
            // no hit, so we just retun 0
            0.0
        }
    }

    // returns the direction from the origin to a random point on this quad
    // TODO: deal with transformations

    fn random(&self, origin: Vec3, rng: &mut SmallRng) -> Vec3 {
        let r0 = rng.gen_range(0.0..1.0);
        let r1 = rng.gen_range(0.0..1.0);

        let p: Vec3 = self.q + (self.u * r0) + (self.v * r1);
        p - origin
    }
}
