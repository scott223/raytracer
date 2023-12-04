use crate::bvh::BVHSplitMethod;
use crate::elements::{Element, HitRecord, Hittable, JSONElement};
use crate::render::{Camera, Interval, Ray};
use crate::{render::camera::JSONCamera, render::Color};

use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Debug, Copy, Clone)]
pub struct Config {
    pub img_width: f64,
    pub img_height: f64,
    pub sample_batch_size: usize,
    pub max_sample_batches: usize,
    pub min_sample_batches: usize,
    pub max_depth: usize,
    pub sky_color: Color,
    pub pixel_radius: f64,
    pub bvh_split_method: Option<BVHSplitMethod>,
}

impl Default for Config {
    fn default() -> Self {
        let r: f64 = 16.0 / 9.0; //aspect ratio
        let w: f64 = 1024.0; //image width
        let h: f64 = w / r; //image heigth, doing the math as double but casting to int as we cannot have a float number of heigth

        let m: usize = 32; //max depth

        let sky_color = Color::new(3.0 / 255.0, 165.0 / 255.0, 252.0 / 255.0);

        let bvh_split_method = BVHSplitMethod::Mid;

        Config {
            img_width: w,
            img_height: h,
            sample_batch_size: 16,
            max_sample_batches: 2,
            min_sample_batches: 1,
            max_depth: m,
            sky_color,
            bvh_split_method: Some(bvh_split_method),
            pixel_radius: 2.0,
        }
    }
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct JSONScene {
    pub camera: JSONCamera,
    pub elements: Vec<JSONElement>,
}

#[derive(Debug)]
pub struct Scene {
    pub camera: Camera,
    pub elements: Vec<Element>,
}

/*
impl Default for Scene {
    fn default() -> Self {
        let lambertian_1 = Material::Lambertian(Lambertian::new(Color::new(0.3, 0.3, 0.3)));
        //       let lambertian_2 = Material::Lambertian(Lambertian::new(Color::new(0.6, 0.7, 0.4)));

        let metal: Material = Material::Metal(Metal::new(Color::new(0.8, 0.7, 0.7), 0.1));
        let glass: Material = Material::Dielectric(Dielectric::new(1.5));

        Scene {
            camera: Camera::defa,
            elements: vec![
                Element::Sphere(Sphere::new(Vec3::new(0.0, -1.0, -25.0), 2.0, metal)),
                Element::Sphere(Sphere::new(Vec3::new(-3.0, -1.25, -23.0), 1.5, glass)),
                Element::Sphere(Sphere::new(Vec3::new(2.0, -2.0, -20.0), 1.0, lambertian_1)),
                //               Element::Plane(Plane::new(Vec3::new(0.0, -2.5, 0.0), Vec3::new(0.0, -1.0, 0.0), lambertian_2)),
            ],
        }
    }
} */

#[allow(dead_code)]
impl Scene {
    // fn trace_hits
    // find the nearest hit for a ray
    pub fn trace(&self, ray: &Ray) -> Option<HitRecord> {
        self.elements
            .iter()
            .filter_map(|e| e.hit(ray, &mut Interval::new(0.001, f64::MAX)))
            .min_by(|i1, i2| i1.t.partial_cmp(&i2.t).unwrap())
    }
}
