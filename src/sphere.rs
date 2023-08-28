use crate::element::HitRecord;
use crate::element::Hittable;
use crate::materials::Material;
use crate::ray::Ray;
use crate::vec3::Vec3;

#[derive(Debug, Clone, Copy)]
pub struct Sphere {
    pub center: Vec3,
    pub radius: f64,
    pub material: Material,
}

impl Sphere {
    // Creating a new Spehere with a center and a radius
    pub fn new(center: Vec3, radius: f64, material: Material) -> Self {
        Sphere {
            center: center,
            radius: radius,
            material: material,
        }
    }
}

// finding the hits for a given ray
// based on: https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-sphere-intersection.html
impl Hittable for Sphere {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let l = ray.origin - self.center;
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
                return None;
            }

            let nearest_root = if first_root < second_root {
                first_root
            } else {
                second_root
            };

            if nearest_root < t_max && nearest_root > t_min {
                let p = ray.at(nearest_root); //we dont need to find the solution with the discrimant, but can just ask the ray where it was at a given t
                let n = p - self.center;
    
                return Some(HitRecord {
                    t: nearest_root,
                    normal: n.normalized(),
                    point: p,
                    material: self.material,
                });
            }
        }
        None // no hits found
    }
}

#[cfg(test)]
mod tests {
    use crate::sphere::Sphere;
    use crate::vec3::Vec3;
    use crate::ray::{Ray};
    use crate::element::{Hittable};
    use crate::materials::Material;
    use crate::materials::Lambertian;
    use crate::color::Color;
    use assert_approx_eq::assert_approx_eq;

    #[test_log::test]
    fn test_create_sphere() {
        let m1: Material = Material::Lambertian(Lambertian::new(Color::new(1.0, 1.0, 1.0)));
        let s = Sphere::new(Vec3::new(1.0, 2.0, -1.0), 2.0, m1);
        assert_approx_eq!(s.center, Vec3::new(1.0, 2.0, -1.0));
        assert_approx_eq!(s.radius, 2.0);
    }

    #[test_log::test]
    fn test_hit_sphere() {
        let m1: Material = Material::Lambertian(Lambertian::new(Color::new(1.0, 1.0, 1.0)));
        let r = Ray::new(Vec3::new(0.0,0.0,0.0),Vec3::new(0.0, 0.0, -1.0));
        let s: Sphere = Sphere::new(Vec3::new(0.0,0.0,-3.0), 1.0, m1);

        if let Some(hit) = s.hit(&r, 0.0, f64::MAX) {
            assert_eq!(hit.t,2.0);
            // assert_eq!(hit.color,s.color);
            assert_eq!(hit.point, Vec3::new(0.0,0.0,-2.0));
        } 
    }
}
