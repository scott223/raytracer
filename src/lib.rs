// Exteral imports
use element::{HitRecord, Hittable};
use image::RgbImage;
use materials::Scatterable;
use std::error::Error;
use std::time::Instant;
use rayon::prelude::*;

// Specfic imports
pub mod config;
mod camera;
mod color;
mod element;
mod materials;
mod plane;
mod ray;
mod sphere;
mod vec3;
use crate::camera::Camera;
use crate::color::Color;
use crate::config::Config;
use crate::ray::Ray;

pub fn render(c: Config) -> Result<(), Box<dyn Error>> {
    // create a new camera object
    let camera: Camera = Camera::new(c.img_width, c.img_width / c.img_height);

    // create a 1-d vector holding all the pixels, and split into bands for parallel rendering
    let mut pixels = vec![Color::new(0.0, 0.0, 0.0); (c.img_width * c.img_height) as usize];
    let bands: Vec<(usize, &mut [Color])> = pixels.chunks_mut(c.img_width as usize).enumerate().collect();

    // use Rayon parallel iterator to iterate over the bands and render per line
    let start = Instant::now();
    bands.into_par_iter().for_each(|(i, band)| {
        render_line(i as i64, band, &camera, &c);
    });

    // Create a new ImgBuf with width: img_width and height: img_y from the config
    let mut img = RgbImage::new(c.img_width as u32, c.img_height as u32);

    // walk through all the pixels, convert to rgb and write to the img pixel as rgb
    for x in 0..(c.img_width) as usize {
        for y in 0..(c.img_height) as usize {
           let pixel = img.get_pixel_mut(x as u32, y as u32);
            *pixel = pixels[(y*c.img_width as usize)+x].to_rgb();
        }
    }

    println!("Frame time: {}ms", start.elapsed().as_millis());
    img.save("renders/render.png").unwrap();
    log::info!("Finished!");
    Ok(())
}

pub fn render_line(y: i64, band: &mut [Color], camera: &Camera, config: &Config) {
    for x in 0..config.img_width as usize {
        let mut color: Color = Color::new(0.0, 0.0, 0.0);

        for _i in 0..config.samples {
            //get multiple rays for anti alliasing, and add the colors
            let ray = camera.get_prime_ray(x as i64, y);
            color += ray_color(&config, &ray, config.max_depth);
        }

        band[x] = color
        .divide_by_samples(config.samples)
        .linear_to_gamma();
    }
}

fn ray_color(c: &Config, r: &Ray, depth: u8) -> Color {
    let hits = trace_hits(&c, &r, 0.001, f64::MAX); //using a very small (but not zero) t_min to avoid shadow acne

    if depth == 0 {
        // we ran out of depth iterations, so we return black
        return Color::new(0.0, 0.0, 0.0);
    }

    match hits {
        Some(hit) => {
            // we hit something
            let scattered = hit.material.scatter(r, &hit);

            match scattered {
                Some((scattered_ray, albedo)) => match scattered_ray {
                    Some(sr) => {
                        let target_color = ray_color(&c, &sr, depth - 1);
                        return Color::new(
                            albedo.r * target_color.r,
                            albedo.g * target_color.g,
                            albedo.b * target_color.b,
                        );
                    }
                    None => albedo,
                },

                None => {
                    // no scattered ray
                    return Color::new(0.0, 0.0, 0.0);
                }
            }
        }
        None => {
            // we did not hit anything, so we paint the sky
            Color::new(3.0 / 255.0, 165.0 / 255.0, 252.0 / 255.0) // return sky
        }
    }
}

// find the nearest hit for a ray, by looping through all the elements on the scene
fn trace_hits(config: &Config, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
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
    use crate::config::Config;
    use crate::render;
    use std::error::Error;

    #[test_log::test]
    fn test_render_full_scene() -> Result<(), Box<dyn Error>> {
        let conf = Config::default();
        render(conf)?;
        Ok(())
    }
}
