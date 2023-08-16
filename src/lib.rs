// Exteral imports
use std::error::Error;
use image::{RgbImage};
use ray::HitRecord;

// Specfic imports
mod vec3;
mod ray;
mod sphere;
pub mod config;
use crate::vec3::Vec3;
use crate::ray::{Ray, Hittable};
use crate::config::Config;

pub fn render(c: &Config) -> Result<(), Box<dyn Error>> {
    // Calculate the vectors across the horizontal and down the vertical viewport edges. this is basically the coordinate of the two corners
    let viewport_u: Vec3 = Vec3::new(c.viewport_u, 0.0, 0.0);
    let viewport_v: Vec3 = Vec3::new(0.0, -c.viewport_v, 0.0);

    // spacing between each pixel
    let pixel_delta_u: Vec3 = viewport_u / (c.img_width as f64);
    let pixel_delta_v: Vec3 = viewport_v / (c.img_height as f64);

    log::info!("pixel_delta_u: {:.5}, pixel_delta_v: {:.5}", pixel_delta_u.x(), pixel_delta_v.y());
    
    // now we position the viewport in front of the camera, taking into account the camera position, the focal length, and the viewport size (we devide those by 2 to put the camera in the middle)
    let viewport_upper_left = c.camera_center - Vec3::new(0.0,0.0,c.focal_length) - viewport_u/2.0 - viewport_v/2.0;
    let pixel00_loc = viewport_upper_left + (pixel_delta_u + pixel_delta_v) * 0.5;
    
    log::info!("viewport_upper_left: {}, pixel00_loc: {}", viewport_upper_left, pixel00_loc);

    // Create a new ImgBuf with width: img_width and height: img_y from the config file
    let mut img = RgbImage::new(c.img_width, c.img_height);
    
    // Start the actual render and creating the rays
    for x in 0..c.img_width {
        for y in 0..c.img_height {
            
            let pixel_loc = pixel00_loc + (pixel_delta_u * x as f64) + (pixel_delta_v * y as f64);
            let ray_direction = pixel_loc - c.camera_center;
            
            let ray = Ray::new(pixel_loc, ray_direction);
            log::trace!("Ray created for pixel ({}.{}): {}", x, y, ray);

            let pixel = img.get_pixel_mut(x, y);
            *pixel = ray_color(&c, &ray);

        }
    }

    img.save("renders/render.png").unwrap();
    log::info!("Finished!");
    Ok(())

}

fn ray_color(c: &Config, r: &Ray) -> image::Rgb<u8> {

    let hits = trace_hits(&c, &r);

    match hits {
        Some(hit) => {
            
            image::Rgb([
                ((hit.normal.x()+1.0)*0.5*255.0) as u8,
                ((hit.normal.y()+1.0)*0.5*255.0) as u8,
                ((hit.normal.z()+1.0)*0.5*255.0) as u8,
            ]) 

        }
        None => {
            image::Rgb([0,0,0])
        }
    }
}

fn trace_hits (
    c: &Config,
    r: &Ray,
) -> Option<HitRecord> {

    let mut hit_list = None;

    //spheres firsts
    for sphere in &c.spheres {
        if let Some(hit) = sphere.hit(&r) {
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