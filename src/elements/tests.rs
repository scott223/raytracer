use std::f64::consts::PI;
use std::fmt::Debug;
use std::fs::File;
use std::io::Read;

use rand::rngs::SmallRng;
use rand::Rng;

use serde::{Deserialize, Serialize};

extern crate wavefront_obj;
use wavefront_obj::obj;

use crate::linalg::{Mat4, Onb, Vec3, Vec4};
use crate::render::Interval;
use crate::render::Ray;
use crate::{bvh::Aabb, materials::*};

#[cfg(test)]
mod tests {
    use crate::render::Color;
    use crate::elements::*;
    use crate::linalg::Vec3;
    use crate::render::Ray;
    use assert_approx_eq::assert_approx_eq;

    #[test_log::test]
    fn test_create_triangle() {
        let m1: Material = Material::Lambertian(Lambertian::new(Color::new(1.0, 1.0, 1.0)));
        let t: Triangle = Triangle::new(
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 2.0, 0.0),
            Vec3::new(2.0, 2.0, 0.0),
            m1,
        );
        assert_approx_eq!(t.v0, Vec3::new(0.0, 0.0, 0.0));
        assert_approx_eq!(t.v1, Vec3::new(0.0, 2.0, 0.0));
        assert_approx_eq!(t.v2, Vec3::new(2.0, 2.0, 0.0));
    }

    #[test_log::test]
    fn test_hit_triangle() {
        let m1: Material = Material::Lambertian(Lambertian::new(Color::new(1.0, 1.0, 1.0)));
        let t: Triangle = Triangle::new(
            Vec3::new(-2.0, -2.0, -5.0),
            Vec3::new(-2.0, 2.0, -5.0),
            Vec3::new(2.0, 2.0, -5.0),
            m1,
        );

        let r: Ray = Ray::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, -1.0));

        match t.hit(&r, &mut Interval::new(0.0, f64::MAX)) {
            Some(_hit) => {
                //assert_eq!(hit.t, 5.0);
                //assert_eq!(hit.front_face, false)
            }
            _ => {
                //panic!("Triangle should be hit")
            }
        }

        //TODO: rewrite test to new MT algo
    }

    #[test_log::test]
    fn test_create_sphere() {
        let m1: Material = Material::Lambertian(Lambertian::new(Color::new(1.0, 1.0, 1.0)));
        let s: Sphere = Sphere::new(Vec3::new(1.0, 2.0, -1.0), 2.0, m1, Some(false));
        assert_approx_eq!(s.center, Vec3::new(1.0, 2.0, -1.0));
        assert_approx_eq!(s.radius, 2.0);
    }

    #[test_log::test]
    fn test_sphere_aabb() {
        let m1: Material = Material::Lambertian(Lambertian::new(Color::new(1.0, 1.0, 1.0)));
        let s: Sphere = Sphere::new(Vec3::new(0.0, 0.0, 0.0), 1.0, m1, Some(false));

        assert_approx_eq!(s.bounding_box().axis(0).interval_min, -1.0);
        assert_approx_eq!(s.bounding_box().axis(0).interval_max, 1.0);

        assert_approx_eq!(s.bounding_box().axis(1).interval_min, -1.0);
        assert_approx_eq!(s.bounding_box().axis(1).interval_max, 1.0);

        assert_approx_eq!(s.bounding_box().axis(2).interval_min, -1.0);
        assert_approx_eq!(s.bounding_box().axis(2).interval_max, 1.0);

        let s2: Sphere = Sphere::new(Vec3::new(2.0, 1.0, -1.0), 1.0, m1, Some(false));

        assert_approx_eq!(s2.bounding_box().axis(0).interval_min, 1.0);
        assert_approx_eq!(s2.bounding_box().axis(0).interval_max, 3.0);

        assert_approx_eq!(s2.bounding_box().axis(1).interval_min, 0.0);
        assert_approx_eq!(s2.bounding_box().axis(1).interval_max, 2.0);

        assert_approx_eq!(s2.bounding_box().axis(2).interval_min, -2.0);
        assert_approx_eq!(s2.bounding_box().axis(2).interval_max, 0.0);
    }

    #[test_log::test]
    fn test_hit_sphere() {
        let m1: Material = Material::Lambertian(Lambertian::new(Color::new(1.0, 1.0, 1.0)));
        let r: Ray = Ray::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, -1.0));
        let s: Sphere = Sphere::new(Vec3::new(0.0, 0.0, -3.0), 1.0, m1, Some(false));

        match s.hit(&r, &mut Interval::new(0.0, f64::MAX)) {
            Some(hit) => {
                assert_eq!(hit.t, 2.0);
                assert_eq!(hit.point, Vec3::new(0.0, 0.0, -2.0));
            }
            _ => {
                panic!("Sphere should be hit")
            }
        }
    }
}
