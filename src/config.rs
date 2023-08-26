use crate::vec3::Vec3;
use crate::sphere::Sphere;
use crate::plane::Plane;
use crate::element::Element;

#[derive(Debug)]
pub struct Config {
    pub ratio: f64,
    pub img_width: f64,
    pub img_height: f64,
    pub samples: u8,
    pub max_depth: u8,
    pub elements: Vec<Element>,
}

impl Default for Config {
    fn default() -> Self {
        let r: f64 = 16.0/9.0; //aspect ratio
        let w: f64 = 800.0; //image width
        let h: f64 = w/r; //image heigth, doing the math as double but casting to int as we cannot have a float number of heigth

        let s: u8 = 1; //samples
        let m: u8 = 32; //max depth

        Config {
            ratio: r,
            img_width: w,
            img_height: h,
            samples: s,
            max_depth: m,
            elements: vec![
                Element::Sphere(Sphere::new(Vec3::new(0.0, -1.0, -25.0), 2.0)),
                Element::Sphere(Sphere::new(Vec3::new(2.0, -2.0, -20.0), 1.0)),
                //Element::Sphere(Sphere::new(Vec3::new(0.0, -61.0, -35.0), 60.0)),                
                Element::Plane(Plane::new(Vec3::new(0.0, -2.5, 0.0), Vec3::new(0.0, -1.0, 0.0))),
            ],
        }
    }
}