use crate::vec3::Vec3;
use crate::sphere::Sphere;

pub struct Config {
    pub ratio: f64,
    pub img_width: f64,
    pub img_height: f64,
    pub samples: u8,
    pub spheres: Vec<Sphere>,
}

impl Default for Config {
    fn default() -> Self {
        let r: f64 = 16.0/9.0; //aspect ratio
        let w: f64 = 800.0; //image width
        let h: f64 = w/r; //image heigth, doing the math as double but casting to int as we cannot have a float number of heigth

        let s: u8 = 4; //samples

        Config {
            ratio: r,
            img_width: w,
            img_height: h,
            samples: s,
            spheres: vec![
                Sphere::new(Vec3::new(0.0, 0.0, -25.0), 2.0),
                Sphere::new(Vec3::new(2.0, -2.0, -18.0), 1.0),
            ],
        }
    }
}