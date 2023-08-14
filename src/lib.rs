// Exteral imports
use std::error::Error;
use image::{Pixel, Rgb, RgbImage};
use std::cmp::PartialEq;
use log::info;

// Specfic imports
mod vec3;
mod ray;
pub mod config;
use crate::vec3::Vec3;
use crate::ray::Ray;
use crate::config::Config;

pub fn render(c: &Config) -> Result<(), Box<dyn Error>> {
    let o: Vec3 = Vec3::new(0 as f64,0 as f64,0 as f64);

    // Calculate the vectors across the horizontal and down the vertical viewport edges. this is basically the coordinate of the two corners
    let viewport_u: Vec3 = Vec3::new(c.viewport_u, 0.0, 0.0);
    let viewport_v: Vec3 = Vec3::new(0.0, -c.viewport_v, 0.0);

    // spacing between each pixel
    let pixel_delta_u: Vec3 = viewport_u / (c.img_width as f64);
    let pixel_delta_v: Vec3 = viewport_v / (c.img_height as f64);

    log::info!("pixel_delta_u: {}, pixel_delta_v: {}", pixel_delta_u.x(), pixel_delta_v.y());
    
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
            
            let ray = Ray::new(pixel_loc, -ray_direction);

            log::trace!("Ray created for pixel ({}.{}): {}", x, y, ray);
            let pixel = img.get_pixel_mut(x, y);

            if (hit_sphere(Vec3::new(0.0,0.0, -3.0), 2.0, &ray)) {
                *pixel = image::Rgb([195,31,31]);
                log::trace!("Ray hitted the sphere!");
            } else {
                *pixel = image::Rgb([0,0,0]);
            }

        }
    }

    img.save("renders/render.png").unwrap();

    log::info!("Finished!");

    Ok(())
}

fn hit_sphere(center: Vec3, radius: f64, ray: &Ray) -> bool {
    let oc = ray.origin - center;
    let a = ray.direction.length_squared();
    let half_b = oc.dot(&ray.direction);
    let c = oc.length_squared() - radius * radius;
    let dscr = (half_b * half_b) - (a * c);

    if (dscr >= 0.0) {
        true
    } else {
        false
    }
}

//bool hit_sphere(const point3& center, double radius, const ray& r) {
//    vec3 oc = r.origin() - center;
//    auto a = dot(r.direction(), r.direction());
//    auto b = 2.0 * dot(oc, r.direction());
//    auto c = dot(oc, oc) - radius*radius;
 //   auto discriminant = b*b - 4*a*c;
 //   return (discriminant >= 0);
//}

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