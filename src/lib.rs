// Rust imports
use image::RgbImage;
use rayon::prelude::*;
use std::error::Error;
use std::time::Instant;

// Raytracer imports
mod camera;
mod color;
pub mod config;
mod element;
mod materials;
mod plane;
mod ray;
mod sphere;
mod vec3;

use camera::Camera;
use color::Color;
use config::Config;
use config::Scene;
use materials::Scatterable;
use ray::Ray;

// fn render
// the main render function that sets up the camera, creates an 1d vector for the pixels, splits it into bands, calls the band render function and writes to an image file
// applies parallel rendering using Rayon
pub fn render(scene: Scene, config: Config) -> Result<(), Box<dyn Error>> {
    // create a new camera object
    let camera: Camera = Camera::new(&config);

    // create a 1-d vector holding all the pixels, and split into bands for parallel rendering
    let mut pixels =
        vec![Color::new(0.0, 0.0, 0.0); (config.img_width * config.img_height) as usize];
    let bands: Vec<(usize, &mut [Color])> = pixels
        .chunks_mut(config.img_width as usize)
        .enumerate()
        .collect();

    log::info!("Starting render!");

    // use Rayon parallel iterator to iterate over the bands and render per line
    let start = Instant::now();
    bands.into_par_iter().for_each(|(i, band)| {
        render_line(i as i64, band, &camera, &scene, &config);
    });

    // Create a new ImgBuf with width: img_width and height: img_y from the config
    let mut img = RgbImage::new(config.img_width as u32, config.img_height as u32);

    // walk through all the pixels and write to the img pixel as rgb
    for x in 0..(config.img_width) as usize {
        for y in 0..(config.img_height) as usize {
            let pixel = img.get_pixel_mut(x as u32, y as u32);
            *pixel = pixels[(y * config.img_width as usize) + x].to_rgb(); // to_rgb applies clamp
        }
    }

    // save the image as a png
    log::info!("Frame time: {}ms", start.elapsed().as_millis());
    img.save("renders/render.png").unwrap();

    // Closing
    log::info!("Finished!");
    Ok(())
} // fn render

// fn render_line
// takes a band of pixels (a single horizontal line along a fixed y coordinate) and renders the pixels on that band
pub fn render_line(y: i64, band: &mut [Color], camera: &Camera, scene: &Scene, config: &Config) {
    for x in 0..config.img_width as usize {
        // start with black color
        let mut color: Color = Color::new(0.0, 0.0, 0.0);

        // loop through all the anti aliasing samples
        for _i in 0..config.samples {
            // get multiple rays for anti alliasing, and add the colors
            let ray = camera.get_prime_ray(x as i64, y);
            color += ray_color(&scene, &config, &ray, config.max_depth);
        }

        // set pixel color, but first divide by the number of samples to get the average, and to the gamma correction
        band[x] = color.divide_by_samples(config.samples).linear_to_gamma();
    }
} // fn render_line

// fn ray_color
// finds the color of the ray, by checking for hits and using the material of that hit + scattered/reflected rays to establish the color
// limit the number of rays it will scatter/reflect by setting depth (e.g. to 32)
fn ray_color(scene: &Scene, config: &Config, ray: &Ray, depth: usize) -> Color {
    if depth == 0 {
        // we ran out of depth iterations, so we return black
        return Color::new(0.0, 0.0, 0.0);
    }

    match scene.trace(&ray) {
        Some(hit) => {
            // we hit something

            match hit.material.scatter(ray, &hit) {
                // lets start with the scatter
                Some((scattered_ray, albedo)) => match scattered_ray {
                    // there is a scattered ray, so lets get the color of that ray
                    Some(sr) => {
                        // call the ray_color function again, but decrease the depth by 1 so that we dont run into an infinite loop
                        let target_color = ray_color(&scene, &config, &sr, depth - 1);

                        // return the color, by applying the albedo to the color of the scattered ray (albedo is here defined the amount of color not absorbed)
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
                    // return black
                    return Color::new(0.0, 0.0, 0.0);
                }
            }
        }
        None => {
            // we did not hit anything, so we paint the sky
            Color::new(3.0 / 255.0, 165.0 / 255.0, 252.0 / 255.0) // return sky
        }
    }
} // fn ray_color

#[cfg(test)]
mod tests {
    use crate::config::Config;
    use crate::config::Scene;
    use crate::render;
    use std::error::Error;

    #[test_log::test]
    fn test_render_full_scene() -> Result<(), Box<dyn Error>> {
        let config: Config = Config::default();
        let scene: Scene = Scene::default();
        render(scene, config)?;
        Ok(())
    }
}