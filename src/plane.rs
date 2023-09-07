use crate::element::HitRecord;
use crate::element::Hittable;
use crate::materials::Material;
use crate::ray::Ray;
use crate::vec3::Vec3;

use serde::{Serialize, Deserialize};

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct Plane {
    pub origin: Vec3,
    pub normal: Vec3,
    pub material: Material,
}

// Creating a new Plane with an origin and a normal (unit vector) (planes have infinite size)
impl Plane {
    pub fn new(origin: Vec3, normal: Vec3, material: Material) -> Self {
        Plane {
            origin: origin,
            normal: normal.normalized(),
            material: material,
        }
    }
}
// finding the hits for a given ray
// return a HitRecord, that contains the position on the ray, the point of the hit & the material
// https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-plane-and-ray-disk-intersection.html
impl Hittable for Plane {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        // get the denominator
        let denom = self.normal.dot(&ray.direction);

        // check if larger than zero, and positive
        if denom > 1e-6 {
            let v = self.origin - ray.origin;
            let distance = v.dot(&self.normal) / denom;

            if distance >= 0.0 {
                if distance > t_min && distance < t_max {
                    let p = ray.at(distance);
                    let hit = HitRecord {
                        t: distance,
                        normal: -self.normal, //we need a minus here to get the defraction working, not sure why.....
                        point: p,
                        material: self.material,
                    };
                    return Some(hit);
                }
            }
        }
        None // no hits found
    }
}

#[cfg(test)]
mod tests {
    use crate::color::Color;
    use crate::element::Hittable;
    use crate::materials::Lambertian;
    use crate::materials::Material;
    use crate::plane::Plane;
    use crate::ray::Ray;
    use crate::vec3::Vec3;
    use assert_approx_eq::assert_approx_eq;

    #[test_log::test]
    fn test_create_plane() {
        let m1: Material = Material::Lambertian(Lambertian::new(Color::new(1.0, 1.0, 1.0)));
        let p = Plane::new(Vec3::new(0.0, -2.5, 0.0), Vec3::new(0.0, -2.0, 0.0), m1);
        assert_approx_eq!(p.origin, Vec3::new(0.0, -2.5, 0.0));
        assert_approx_eq!(p.normal, Vec3::new(0.0, -1.0, 0.0));
    }

    #[test_log::test]
    fn test_hit_plane() {
        let m1: Material = Material::Lambertian(Lambertian::new(Color::new(1.0, 1.0, 1.0)));
        let r = Ray::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0));
        let p = Plane::new(Vec3::new(0.0, 0.0, 5.0), Vec3::new(0.0, 0.0, 1.0), m1);

        if let Some(hit) = p.hit(&r, 0.0, f64::MAX) {
            assert_eq!(hit.t, 5.0);
            assert_eq!(hit.point, Vec3::new(0.0, 0.0, 5.0));
        }

        let r = Ray::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, -1.0, 0.0));
        let p = Plane::new(Vec3::new(0.0, -2.5, 0.0), Vec3::new(0.0, -1.0, 0.0), m1);

        if let Some(hit) = p.hit(&r, 0.0, f64::MAX) {
            assert_eq!(hit.t, 2.5);
            assert_eq!(hit.point, Vec3::new(0.0, -2.5, 0.0));
        }
    }
}
