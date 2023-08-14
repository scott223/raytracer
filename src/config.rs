use crate::vec3::Vec3;
use crate::sphere::Sphere;

use log::info;

pub struct Config {
    pub ratio: f64,
    pub img_width: u32,
    pub img_height: u32,
    pub viewport_u: f64,
    pub viewport_v: f64,
    pub focal_length: f64,
    pub camera_center: Vec3,
    pub spheres: Vec<Sphere>,
}

impl Default for Config {
    fn default() -> Self {
        let r: f64 = 16.0/9.0; //aspect ratio
        let w: f64 = 400.0; //image width
        let h: u32 = (w/r) as u32; //image heigth, doing the math as double but casting to int as we cannot have a float number of heigth

        let real_r: f64 = w/h as f64; // getting the real ration, which might be slightly off the desired one
 
        let vv: f64 = 2.0; // viewport height 
        let vu: f64 = vv * real_r; // calculating viewport width

        log::info!("Desired aspect ration: {:.3}, actual aspect ratio: {:.3}. img_width: {}, img_heigth:{}, v_u (viewport width): {}, v_v (viewport heigth): {}", r, real_r, w, h, vu, vv);

        Config {
            ratio: r,
            img_width: w as u32,
            img_height: h as u32,
            viewport_u: vu as f64,
            viewport_v: vv as f64,
            focal_length: 1.0,
            camera_center: Vec3::new(0.0, 0.0, 0.0),
            spheres: vec![
                Sphere::new(Vec3::new(0.0, 0.0, -5.0),2.0),
                Sphere::new(Vec3::new(2.0, -3.0, -5.0),1.0),
            ],
        }
    }
}