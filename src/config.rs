use crate::color::Color;
use crate::ray::Ray;
use crate::vec3::Vec3;
use crate::elements::*;
use crate::materials::*;

use serde::{Serialize, Deserialize};

#[derive(Serialize, Deserialize, Debug)]
pub struct Config {
    pub img_width: f64,
    pub img_height: f64,
    pub camera_center: Vec3,
    pub camera_look_at: Vec3,
    pub camera_fov_vertical: f64,
    pub samples: usize,
    pub max_depth: usize,
    pub sky_color: Color,
}

impl Default for Config {
    fn default() -> Self {
        let r: f64 = 16.0/9.0; //aspect ratio
        let w: f64 = 1024.0; //image width
        let h: f64 = w/r; //image heigth, doing the math as double but casting to int as we cannot have a float number of heigth
        let fov: f64 = 90.0;

        let s: usize = 1; //samples
        let m: usize = 32; //max depth

        let sky_color = Color::new(3.0 / 255.0, 165.0 / 255.0, 252.0 / 255.0);

        Config {
            img_width: w,
            img_height: h,
            camera_center: Vec3::new(0.0, 0.0, 0.0),
            camera_look_at: Vec3::new(0.0, 0.0, -3.0),
            camera_fov_vertical: fov,
            samples: s,
            max_depth: m,
            sky_color,
        }
    }
}

#[derive(Serialize, Deserialize, Debug)]
pub struct Scene {
    pub elements: Vec<Element>,

}

impl Default for Scene {
    fn default() -> Self {

        let lambertian_1 = Material::Lambertian(Lambertian::new(Color::new(0.3, 0.3, 0.3)));
        let lambertian_2 = Material::Lambertian(Lambertian::new(Color::new(0.6, 0.7, 0.4)));

        let metal: Material = Material::Metal(Metal::new(Color::new(0.8, 0.7, 0.7), 0.1));
        let glass: Material = Material::Dielectric(Dielectric::new(1.5));

        Scene {
            elements: vec![
                Element::Sphere(Sphere::new(Vec3::new(0.0, -1.0, -25.0), 2.0, metal)),
                Element::Sphere(Sphere::new(Vec3::new(-3.0, -1.25, -23.0), 1.5, glass)),
                Element::Sphere(Sphere::new(Vec3::new(2.0, -2.0, -20.0), 1.0, lambertian_1)),             
                Element::Plane(Plane::new(Vec3::new(0.0, -2.5, 0.0), Vec3::new(0.0, -1.0, 0.0), lambertian_2)),
            ],
        }
    }
}

impl Scene {

    // fn trace_hits
    // find the nearest hit for a ray
    pub fn trace(&self, ray: &Ray) -> Option<HitRecord> {
        self.elements
            .iter()
            .filter_map(|e| e.hit(&ray, 0.001, f64::MAX))
            .min_by(|i1, i2 | i1.t.partial_cmp(&i2.t).unwrap())

    }
}