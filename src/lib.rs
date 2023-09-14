// Rust imports
use image::{ImageBuffer, Rgb, RgbImage};
use rayon::prelude::*;
use std::error::Error;
use std::time::Instant;

// Raytracer imports
pub mod config;

mod camera;
mod color;
mod elements;
mod materials;
mod ray;
mod vec3;

use camera::Camera;
use color::Color;
use config::{Config, Scene};
use materials::Scatterable;
use ray::Ray;

// fn render
// the main render function that sets up the camera, creates an 1d vector for the pixels, splits it into bands, calls the band render function and writes to an image file
// applies parallel rendering using Rayon
pub fn render(
    scene: Scene,
    config: Config,
) -> Result<ImageBuffer<Rgb<u8>, Vec<u8>>, Box<dyn Error>> {
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

    // create a new ImgBuf with width: img_width and height: img_y from the config
    let mut img = RgbImage::new(config.img_width as u32, config.img_height as u32);

    // walk through all the pixels and write to the img pixel as rgb
    for x in 0..(config.img_width) as usize {
        for y in 0..(config.img_height) as usize {
            let pixel = img.get_pixel_mut(x as u32, y as u32);
            *pixel = pixels[(y * config.img_width as usize) + x]
                .linear_to_gamma()
                .to_rgb(); // to_rgb applies clamp
        }
    }

    // save the image as a png
    log::info!("Render finished in {}ms", start.elapsed().as_millis());
    Ok(img)
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

    // we trace the ray to see what it will hit. the scene.trace iterates through all the objects in the scene and calls the hittable->hit function to determine if a hit is made
    match scene.trace(&ray) {
        Some(hit) => {
            // we hit something
            // we start with scattered rays (we assume every object has scattered rays, although in some materials (like metal) its actually a reflected or refracted (glass) ray)
            match hit.material.scatter(ray, &hit) {

                Some((scattered_ray, albedo)) => {
                    // see if there is a ray returned
                    match scattered_ray {
                        Some(sr) => {
                            // there is a scattered ray, so lets get the color of that ray
                            // call the ray_color function again, now with the reflected ray but decrease the depth by 1 so that we dont run into an infinite loop
                            let target_color = ray_color(&scene, &config, &sr, depth - 1);

                            // return the color, by applying the albedo to the color of the scattered ray (albedo is here defined the amount of color not absorbed)
                            return Color::new(
                                albedo.r * target_color.r,
                                albedo.g * target_color.g,
                                albedo.b * target_color.b,
                            );
                        }
                        None => {
                            // there is no ray, so just return the albedo of the hittable object (this could be a light?)
                            albedo
                        }
                    }
                }
                None => {
                    // no scattered ray
                    // return black
                    return Color::new(0.0, 0.0, 0.0);
                }
            }
        }
        None => {
            // we did not hit anything, so we return the color of the sky
            config.sky_color
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
