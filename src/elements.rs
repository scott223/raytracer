use std::fmt::Debug;
use crate::{materials::*, aabb::Aabb};
use crate::ray::Ray;
use crate::vec3::Vec3;
use crate::interval::Interval;

use serde::{Deserialize, Serialize};

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
}

// enum for all the different elements
#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub enum Element {
    Sphere(Sphere),
    Quad(Quad),
//    Plane(Plane),
}

// matching the hit function with the hittable trait for each type of element
impl Hittable for Element {
    fn hit(&self, ray: &Ray, ray_t: &mut Interval) -> Option<HitRecord> {
        match *self {
            Element::Sphere(ref s) => s.hit(ray, ray_t),
            Element::Quad(ref q) => q.hit(ray, ray_t),
 //           Element::Plane(ref p) => p.hit(ray, t_min, t_max),
        }
    }

    fn bounding_box(&self) -> Aabb {
        match *self {
            Element::Sphere(ref s) => s.bounding_box(),
            Element::Quad(ref q) => q.bounding_box(),
        }
    }
}

// Plane element
#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct Plane {
    pub origin: Vec3,
    pub normal: Vec3,
    pub material: Material,
}

// Creating a new Plane with an origin and a normal (unit vector) (planes have infinite size)
impl Plane {
    #[allow(dead_code)]
    pub fn new(origin: Vec3, normal: Vec3, material: Material) -> Self {
        Plane {
            origin,
            normal: normal.normalized(),
            material,
        }
    }
}

// finding the hits for a given ray
// return a HitRecord, that contains the position on the ray, the point of the hit & the material
// https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-plane-and-ray-disk-intersection.html
impl Hittable for Plane {
    fn hit(&self, ray: &Ray, ray_t: &mut Interval) -> Option<HitRecord> {
        // get the denominator
        let denom = self.normal.dot(&ray.direction);

        // check if larger than zero, and positive
        if denom > 1e-6 {
            let v = self.origin - ray.origin;
            let distance = v.dot(&self.normal) / denom;

            if distance >= 0.0 && distance > ray_t.interval_min && distance < ray_t.interval_max {
                let p = ray.at(distance);
                let hit = HitRecord {
                    t: distance,
                    normal: -self.normal, //we need a minus here to get the defraction working, not sure why.....
                    point: p,
                    front_face: true,
                    material: self.material,
                };
                return Some(hit);
            }
        }
        None // no hits found
    }

    fn bounding_box(&self) -> Aabb {
        // TODO, will need to figure out how to add a bounding box for an infinite plane
        todo!();
    }
}

//Quad element
#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct Quad {
    pub Q: Vec3,
    pub u: Vec3,
    pub v: Vec3,
    pub material: Material,
    // pub bbox: Aabb,
    // TODO initialize the bounding box at object creation!
}



// Creating a new Sphere with a center and a radius
impl Quad {
    pub fn new(Q: Vec3, u: Vec3, v:Vec3, material: Material) -> Self {
        Quad {
            Q,
            u,
            v,
            material,
        }
    }

    //given the hit point in plane coordinates, return none if it is outside the primitive
    pub fn is_interior(a: f64, b: f64) -> Option<Vec<f64>> {
        if (a < 0.0) || (1.0 < a) || (b < 0.0) || (1.0 < b) {
            return None
        }

        // we have this work around as we dont determine the u and v in the hitrecord yet
        let mut result: Vec<f64> = Vec::new();
        result.push(a);
        result.push(b);

        return Some(result);
    }
}

impl Hittable for Quad {
    fn hit(&self, ray: &Ray, ray_t: &mut Interval) -> Option<HitRecord> {
        
        //TODO: we should preload the values vor n, normal and dot in the constructor
        let n: Vec3 = self.u.cross(&self.v);
        let normal: Vec3 = n.normalized();
        let D: f64 = normal.dot(&self.Q);
        let w: Vec3 = n / n.dot(&n);

        let denom = normal.dot(&ray.direction);

        // no hits, as the ray is parallel to the plane
        if denom.abs() < 0.000000001 {
            return None
        }

        // return false if the hit point paramater t is outside the ray interval
        let t: f64 = (D - normal.dot(&ray.origin)) / denom;
        if !ray_t.contains(t) {
            return None;
        }

        let intersection = ray.at(t);
        let planar_hitpoint_vector: Vec3 = intersection - self.Q;
        let alpha: f64 = w.dot(&planar_hitpoint_vector.cross(&self.v));
        let beta: f64 = w.dot(&self.u.cross(&planar_hitpoint_vector));

        match Quad::is_interior(alpha, beta) {
            Some(_result) => {
                // nothing needed here i think
            },
            _ => { 
                //is not in the interior, so break
                return None 
            }
        }

        // we have a hit, so we return a hit record
        return Some(HitRecord {
            t,
            normal, //TODO check normal (should be outward normal)
            point: intersection,
            front_face: ray.direction.dot(&normal) < 0.0,
            material: self.material,
        });

    }
    
    fn bounding_box(&self) -> Aabb {
        Aabb::new_from_points(self.Q, self.Q + self.u + self.v).pad()
    }
}

// Sphere element
#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct Sphere {
    pub center: Vec3,
    pub radius: f64,
    pub material: Material,
    // TODO initialize the bounding box at object creation!
}

// Creating a new Sphere with a center and a radius
impl Sphere {
    pub fn new(center: Vec3, radius: f64, material: Material) -> Self {
        Sphere {
            center,
            radius,
            material,
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
                    normal: if front_face { outward_normal } else { -outward_normal },
                    point: p,
                    front_face: front_face,
                    material: self.material,
                });
            }
        }
        None // no hits found
    }

    // construct an axis aligned bounding box Aabb for a sphere
    fn bounding_box(&self) -> Aabb {
        let rvec = Vec3::new(self.radius, self.radius, self.radius);
        Aabb::new_from_points(self.center - rvec, self.center + rvec)
    }
}

#[cfg(test)]
mod tests {
    use crate::color::Color;
    use crate::elements::*;
    use crate::ray::Ray;
    use crate::vec3::Vec3;
    use assert_approx_eq::assert_approx_eq;

    #[test_log::test]
    fn test_create_plane() {
        let m1: Material = Material::Lambertian(Lambertian::new(Color::new(1.0, 1.0, 1.0)));
        let p: Plane = Plane::new(Vec3::new(0.0, -2.5, 0.0), Vec3::new(0.0, -2.0, 0.0), m1);
        assert_approx_eq!(p.origin, Vec3::new(0.0, -2.5, 0.0));
        assert_approx_eq!(p.normal, Vec3::new(0.0, -1.0, 0.0));
    }

    #[test_log::test]
    fn test_hit_plane() {
        let m1: Material = Material::Lambertian(Lambertian::new(Color::new(1.0, 1.0, 1.0)));
        let r: Ray = Ray::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0));
        let p: Plane = Plane::new(Vec3::new(0.0, 0.0, 5.0), Vec3::new(0.0, 0.0, 1.0), m1);

        if let Some(hit) = p.hit(&r, &mut Interval::new(0.0, f64::MAX)) {
            assert_eq!(hit.t, 5.0);
            assert_eq!(hit.point, Vec3::new(0.0, 0.0, 5.0));
        }

        let r: Ray = Ray::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, -1.0, 0.0));
        let p: Plane = Plane::new(Vec3::new(0.0, -2.5, 0.0), Vec3::new(0.0, -1.0, 0.0), m1);

        if let Some(hit) = p.hit(&r, &mut Interval::new(0.0, f64::MAX)) {
            assert_eq!(hit.t, 2.5);
            assert_eq!(hit.point, Vec3::new(0.0, -2.5, 0.0));
        }
    }

    #[test_log::test]
    fn test_create_sphere() {
        let m1: Material = Material::Lambertian(Lambertian::new(Color::new(1.0, 1.0, 1.0)));
        let s: Sphere = Sphere::new(Vec3::new(1.0, 2.0, -1.0), 2.0, m1);
        assert_approx_eq!(s.center, Vec3::new(1.0, 2.0, -1.0));
        assert_approx_eq!(s.radius, 2.0);
    }

    #[test_log::test]
    fn test_sphere_aabb() {
        let m1: Material = Material::Lambertian(Lambertian::new(Color::new(1.0, 1.0, 1.0)));
        let s: Sphere = Sphere::new(Vec3::new(0.0, 0.0, 0.0), 1.0, m1);

        assert_approx_eq!(s.bounding_box().axis(0).interval_min,-1.0);
        assert_approx_eq!(s.bounding_box().axis(0).interval_max,1.0);

        assert_approx_eq!(s.bounding_box().axis(1).interval_min,-1.0);
        assert_approx_eq!(s.bounding_box().axis(1).interval_max,1.0);

        assert_approx_eq!(s.bounding_box().axis(2).interval_min,-1.0);
        assert_approx_eq!(s.bounding_box().axis(2).interval_max,1.0);

        let s2: Sphere = Sphere::new(Vec3::new(2.0, 1.0, -1.0), 1.0, m1);

        assert_approx_eq!(s2.bounding_box().axis(0).interval_min,1.0);
        assert_approx_eq!(s2.bounding_box().axis(0).interval_max,3.0);

        assert_approx_eq!(s2.bounding_box().axis(1).interval_min,0.0);
        assert_approx_eq!(s2.bounding_box().axis(1).interval_max,2.0);
         
        assert_approx_eq!(s2.bounding_box().axis(2).interval_min,-2.0);
        assert_approx_eq!(s2.bounding_box().axis(2).interval_max,0.0);                   
    }

    #[test_log::test]
    fn test_hit_sphere() {
        let m1: Material = Material::Lambertian(Lambertian::new(Color::new(1.0, 1.0, 1.0)));
        let r: Ray = Ray::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, -1.0));
        let s: Sphere = Sphere::new(Vec3::new(0.0, 0.0, -3.0), 1.0, m1);

        if let Some(hit) = s.hit(&r, &mut Interval::new(0.0, f64::MAX)) {
            assert_eq!(hit.t, 2.0);
            assert_eq!(hit.point, Vec3::new(0.0, 0.0, -2.0));
        }
    }
}
