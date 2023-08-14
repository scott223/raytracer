use crate::vec3::Vec3;
use crate::ray::Ray;
use crate::ray::HitRecord;
use crate::ray::Hittable;
use std::fmt;

#[derive(Debug, Clone, Copy)]
pub struct Sphere {
    pub center: Vec3,
    pub radius: f64,
    pub color: image::Rgb<u8>,
}

impl Sphere {
    // Creating a new Spehere with a center and a radius
    pub fn new(center_val: Vec3, r_val: f64) -> Self {
        Sphere{
            center: center_val,
            radius: r_val,
            color: image::Rgb([80,31,31]),
        }
    }
}

impl Hittable for Sphere {

    // finding the hits for a given ray
    // based on: https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-sphere-intersection.html

    fn hit(&self, ray: &Ray) -> Option<HitRecord> {
        let l = ray.origin - self.center;
        let a = ray.direction.dot(&ray.direction);
        let b = 2.0 * ray.direction.dot(&l);
        let c = l.dot(&l) - (self.radius * self.radius);

        let discr: f64 = b * b - 4.0 * a * c;
        
        if (discr >= 0.0) { // there is one or two solutions
            let first_root = -0.5 * (b + discr.sqrt());
            let second_root = -0.5 * (b - discr.sqrt()); 
            
            for root in [first_root, second_root].iter() {
                let p = ray.at(*root); //we dont need to find the solution with the discrimant, but can just ask the ray where it was at a given t
                return Some(HitRecord { 
                    t: *root, 
                    point: p, 
                    color: self.color 
                })
            }
        }
        None // no hits found
    }
}

#[cfg(test)]
mod tests {
    use assert_approx_eq::assert_approx_eq;
    use crate::vec3::Vec3;
    use crate::sphere::Sphere;

    #[test_log::test]
    fn test_create_sphere() {
        let s = Sphere::new(Vec3::new(1.0, 2.0, -1.0), 2.0);
        assert_approx_eq!(s.center, Vec3::new(1.0, 2.0, -1.0));
        assert_approx_eq!(s.radius, 2.0);
    } 
}