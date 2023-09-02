use crate::color::Color;
use crate::vec3::Vec3;
use crate::sphere::Sphere;
use crate::plane::Plane;
use crate::element::Element;
use crate::materials::Material;
use crate::materials::Lambertian;
use crate::materials::Metal;

#[derive(Debug)]
pub struct Config {
    pub ratio: f64,
    pub img_width: f64,
    pub img_height: f64,
    pub samples: usize,
    pub max_depth: usize,
    pub elements: Vec<Element>,
}

impl Default for Config {
    fn default() -> Self {
        let r: f64 = 16.0/9.0; //aspect ratio
        let w: f64 = 1024.0; //image width
        let h: f64 = w/r; //image heigth, doing the math as double but casting to int as we cannot have a float number of heigth

        let s: usize = 64; //samples
        let m: usize = 32; //max depth

        let m1 = Material::Lambertian(Lambertian::new(Color::new(0.3, 0.3, 0.3)));
        let m2 = Material::Lambertian(Lambertian::new(Color::new(0.6, 0.7, 0.4)));

        let m3: Material = Material::Metal(Metal::new(Color::new(0.8, 0.7, 0.7), 0.1));

        Config {
            ratio: r,
            img_width: w,
            img_height: h,
            samples: s,
            max_depth: m,
            elements: vec![
                Element::Sphere(Sphere::new(Vec3::new(0.0, -1.0, -25.0), 2.0, m3)),
                Element::Sphere(Sphere::new(Vec3::new(2.0, -2.0, -20.0), 1.0, m1)),             
                Element::Plane(Plane::new(Vec3::new(0.0, -2.5, 0.0), Vec3::new(0.0, -1.0, 0.0), m2)),
            ],
        }
    }
}