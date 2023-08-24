// Exteral imports
use std::error::Error;
use image::{RgbImage, Pixel};
use element::{Hittable, HitRecord};

// Specfic imports
mod vec3;
mod ray;
mod sphere;
mod plane;
mod element;
mod camera;
pub mod config;
use crate::vec3::Vec3;
use crate::ray::Ray;
use crate::config::Config;
use crate::camera::Camera;

pub fn render(c: &Config) -> Result<(), Box<dyn Error>> {

    let camera: Camera = Camera::new(c.img_width, c.img_width/c.img_height);

    // Create a new ImgBuf with width: img_width and height: img_y from the config file
    let mut img = RgbImage::new(c.img_width as u32, c.img_height as u32);
    
    // Start the actual render and creating the rays
    for x in 0..c.img_width as i64 {
        for y in 0..c.img_height as i64 {
            
            let pixel = img.get_pixel_mut(x as u32, y as u32);

            let mut color = Vec3::new(0.0,0.0,0.0);

            for i in 0..c.samples {
                let ray = camera.get_prime_ray(x, y);
                color = add_clamp(color, ray_color(&c, &ray, c.max_depth), c.samples);
            }

            *pixel = image::Rgb([color.x() as u8, color.y() as u8, color.z() as u8]);

        }
    }

    img.save("renders/render.png").unwrap();
    log::info!("Finished!");
    Ok(())

}

fn add_clamp(color: Vec3, color_to_add: Vec3, samples: u8) -> Vec3 {
    let mut r: f64 = color.x() + color_to_add.x()/((samples) as f64);
    let mut g: f64 = color.y() + color_to_add.y()/((samples) as f64);    
    let mut b: f64 = color.z() + color_to_add.z()/((samples) as f64);

    if r > 255.0 { r = 255.0; }  
    if g > 255.0 { g = 255.0; }
    if b > 255.0 { b = 255.0; }

    Vec3::new(r,g,b)
}

fn ray_color(c: &Config, r: &Ray, depth: u8) -> Vec3 {

    let hits = trace_hits(&c, &r, 0.001, f64::MAX); //using a very small (but not zero) t_min to avoid shadow acne

    if depth == 0 {
        return Vec3::new(0.0,0.0,0.0);
    }

    match hits {
        Some(hit) => {

            let new_direction = hit.normal + Vec3::new_random_unit_vector(); //lambertian distribution
            let color = ray_color(&c, &Ray::new(hit.point,new_direction), depth-1) * 0.5;
            //log::info!("color: {:?}", color);
            return color;
            
            //let color = Vec3::new(
            //    ((hit.normal.x()+1.0)*0.5*255.0) as f64,
            //    ((hit.normal.y()+1.0)*0.5*255.0) as f64,
            //    ((hit.normal.z()+1.0)*0.5*255.0) as f64,
            //);

            //color

        }
        None => {
            Vec3::new(179.0,236.0,255.0)
        }
    }
}

fn trace_hits (
    config: &Config,
    ray: &Ray,
    t_min: f64,
    t_max: f64,
) -> Option<HitRecord> {

    let mut hit_list = None;

    let mut closest_so_far: f64 = t_max;

    for element in &config.elements {
        if let Some(hit) = element.hit(&ray, t_min, closest_so_far) {
            closest_so_far = hit.t;
            hit_list = Some(hit);
        } 
    }

    hit_list
}

#[cfg(test)]
mod tests {
    use std::error::Error;
    use crate::render;
    use crate::config::Config;
    
    #[test_log::test]
    fn test_render_full_scene() -> Result <(), Box<dyn Error>> {
        let conf = Config::default(); 
        render(&conf)?;
        Ok(())
    }
}