use std::f64::consts::PI;

use rand::{Rng, rngs::SmallRng};

use serde::Deserialize;
use serde::Serialize;

use crate::bvh::Aabb;
use crate::linalg::{Mat4, Onb, Vec3, Vec4};
use crate::materials::Material;
use crate::render::{Ray, Interval};

use super::Element;
use super::HitRecord;
use super::Hittable;
use super::Transpose;

// Sphere element
// simplified JSON representation of the Sphere
#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct JSONSphere {
    pub center: Vec3,
    pub radius: f64,
    pub material: Material,
    pub attractor: Option<bool>,
    pub transpose: Option<Transpose>,
}

impl JSONSphere {
    pub fn add_as_element(&self, objects: &mut Vec<Element>) {
        let mut center: Vec3 = self.center;

        // Transpose transformation
        if let Some(t) = self.transpose {
            let tm = Mat4::transpose(t.x, t.y, t.z);
            center = (tm * Vec4::new_from_vec3(center, 1.0)).to_vec3();
        }

        let object = Element::Sphere(Sphere::new(
            center,
            self.radius,
            self.material,
            self.attractor,
        ));
        objects.push(object);
    }
}

// the actual sphere element
#[derive(Debug, Clone, Copy)]
pub struct Sphere {
    pub center: Vec3,
    pub radius: f64,
    pub material: Material,
    pub bbox: Aabb,
    pub attractor: Option<bool>,
}

// Creating a new Sphere with a center and a radius
impl Sphere {
    // create a new sphere based on the JSON object
    pub fn new_from_json_object(json_sphere: JSONSphere) -> Self {
        Sphere::new(
            json_sphere.center,
            json_sphere.radius,
            json_sphere.material,
            json_sphere.attractor,
        )
    }

    // create a new sphere
    pub fn new(center: Vec3, radius: f64, material: Material, attractor: Option<bool>) -> Self {
        let rvec = Vec3::new(radius, radius, radius);
        let bbox = Aabb::new_from_points(center - rvec, center + rvec);

        Sphere {
            center,
            radius,
            material,
            bbox,
            attractor,
        }
    }

    pub fn is_attractor(&self) -> bool {
        if let Some(true) = self.attractor {
            true
        } else {
            false
        }
    }
}

// implementing hte hittable traits for the Sphere
impl Hittable for Sphere {
    // finding the hits for a given ray
    // based on: https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-sphere-intersection.html
    fn hit(&self, ray: &Ray, ray_t: &mut Interval) -> Option<HitRecord> {
        let l = ray.origin - self.center;

        // solve the quadratic equation
        let a = ray.direction.dot(&ray.direction);
        let b = 2.0 * ray.direction.dot(&l);
        let c = l.dot(&l) - (self.radius * self.radius);

        let discr: f64 = b * b - 4.0 * a * c;

        if discr >= 0.0 {
            // there is one or two solutions

            let q: f64 = if b > 0.0 {
                -0.5 * (b + discr.sqrt())
            } else {
                -0.5 * (b - discr.sqrt())
            };
            let first_root: f64 = q / a;
            let second_root: f64 = c / q;

            if first_root < 0.0 && second_root < 0.0 {
                // only negative solutions
                return None;
            }

            let nearest_root = if first_root < second_root {
                first_root
            } else {
                second_root
            };

            if nearest_root < ray_t.interval_max && nearest_root > ray_t.interval_min {
                let p = ray.at(nearest_root); //we dont need to find the solution with the discrimant, but can just ask the ray where it was at a given t
                let outward_normal = ((p - self.center) / self.radius).normalized();
                let front_face = ray.direction.dot(&outward_normal) < 0.0;

                return Some(HitRecord {
                    t: nearest_root,
                    normal: if front_face {
                        outward_normal
                    } else {
                        -outward_normal
                    },
                    point: p,
                    front_face,
                    material: self.material,
                });
            }
        }
        None // no hits found
    }

    // return the axis aligned bounding box Aabb for a sphere
    fn bounding_box(&self) -> Aabb {
        self.bbox
    }

    // returns the probability distribution factor for a given ray
    // source: ray tracing in one weekend
    fn pdf_value(&self, origin: Vec3, direction: Vec3) -> f64 {
        // This method only works for stationary spheres

        let ray = Ray::new(origin, direction);
        if let Some(_hit) = self.hit(&ray, &mut Interval::new(0.001, f64::INFINITY)) {
            let cos_theta_max =
                (1. - self.radius * self.radius / (self.center - origin).length_squared()).sqrt();
            let solid_angle = 2.0 * PI * (1.0 - cos_theta_max);

            //println!("did hit");
            1. / solid_angle
        } else {
            //no hit, return 0 as pdf
            0.
        }
    }

    // returns a random point on the sphere
    // source: raytracing in one weekend
    fn random(&self, origin: Vec3, rng: &mut SmallRng) -> Vec3 {
        let direction: Vec3 = self.center - origin;
        let uvw: Onb = Onb::build_from_w(direction);
        uvw.local_vec(random_to_sphere(
            self.radius,
            direction.length_squared(),
            rng,
        ))
    }
}

pub fn random_to_sphere(radius: f64, distance_squared: f64, rng: &mut SmallRng) -> Vec3 {
    let r1: f64 = rng.gen_range(0.0..1.0);
    let r2: f64 = rng.gen_range(0.0..1.0);
    let z: f64 = 1.0 + r2 * ((1.0 - radius * radius / distance_squared).sqrt() - 1.0);

    let phi: f64 = 2.0 * PI * r1;
    let x: f64 = phi.cos() * (1.0 - z * z).sqrt();
    let y: f64 = phi.sin() * (1.0 - z * z).sqrt();

    return Vec3::new(x, y, z);
}
