// Exteral imports
use std::error::Error;
use image::RgbImage;
use element::{Hittable, HitRecord};

// Specfic imports
mod vec3;
mod ray;
mod sphere;
mod plane;
mod color;
mod element;
mod camera;
pub mod config;
use crate::vec3::Vec3;
use crate::ray::Ray;
use crate::color::Color;
use crate::config::Config;
use crate::camera::Camera;

pub fn render(c: &Config) -> Result<(), Box<dyn Error>> {
    let camera: Camera = Camera::new(c.img_width, c.img_width/c.img_height);

    // Create a new ImgBuf with width: img_width and height: img_y from the config
    let mut img = RgbImage::new(c.img_width as u32, c.img_height as u32);
    
    // Start the actual render and creating the rays, by looping through all the image pixels
    for x in 0..c.img_width as i64 {
        for y in 0..c.img_height as i64 {
            
            let pixel = img.get_pixel_mut(x as u32, y as u32);
            let mut color: Color = Color::new(0.0, 0.0, 0.0);

            for _i in 0..c.samples { //get multiple rays for anti alliasing, and add the colors
                let ray = camera.get_prime_ray(x, y);
                color += ray_color(&c, &ray, c.max_depth);
            }

            *pixel = color.divide_by_samples(c.samples).linear_to_gamma().to_rgb(); //to_rgb function also clamps to max 255 to avoid overflow
        }
    }

    img.save("renders/render.png").unwrap();
    log::info!("Finished!");
    Ok(())
}

fn ray_color(c: &Config, r: &Ray, depth: u8) -> Color {
    let hits = trace_hits(&c, &r, 0.001, f64::MAX); //using a very small (but not zero) t_min to avoid shadow acne

    if depth == 0 { // we ran out of depth iterations, so we return black
        return Color::new(0.0,0.0,0.0);
    }

    match hits {
        Some(hit) => { // we hit something, so cast a new ray
            let new_direction = hit.normal + Vec3::new_random_unit_vector(); //lambertian distribution
            ray_color(&c, &Ray::new(hit.point,new_direction), depth-1) * 0.3
        }
        None => {
            Color::new(3.0/255.0,165.0/255.0,252.0/255.0) // return sky
        }
    }
}

// find the nearest hit for a ray, by looping through all the elements on the scene
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