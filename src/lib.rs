// Rust imports
use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};

use indicatif::ProgressBar;

use rayon::prelude::*;
use std::error::Error;
use std::f64::consts::PI;
use std::time::Instant;

// Raytracer imports
pub mod config;

mod aabb;
mod bhv;
mod camera;
pub mod color;
pub mod elements;
pub mod interval;
pub mod materials;
pub mod ray;
pub mod vec3;
mod mat4;
mod onb;

use bhv::BHVNode;
use camera::Camera;
use color::Color;
use config::{Config, JSONScene, Scene};
use elements::{Element, Hittable};
use interval::Interval;
use materials::{Emmits, Scatterable};
use ray::Ray;

use crate::elements::{JSONElement, Triangle};
use crate::elements::{Sphere, Quad};

use crate::materials::{Material, Lambertian};
use crate::vec3::Vec3;

// fn render
// the main render function that sets up the camera, creates an 1d vector for the pixels, splits it into bands, calls the band render function and writes to an image file
// applies parallel rendering using Rayon
pub fn render(
    json_scene: JSONScene,
    config: Config,
) -> Result<Vec<Color>, Box<dyn Error>> {
    // create a new camera object
    let camera: Camera = Camera::new(&config);

    // we loop through the JSON Scene, and create a new actual object for each simplified JSON object
    // this allows us to pre-calculate various things like the bbox or commonly used numbers in the cconstructor
    let mut objects: Vec<Element> = Vec::new();

    for json_element in json_scene.elements {
        match json_element {
            JSONElement::JSONQuad(q) => q.add_as_element(&mut objects),
            JSONElement::JSONSphere(s) => s.add_as_element(&mut objects),
            JSONElement::JSONTriangle(t) => t.add_as_element(&mut objects),
            JSONElement::JSONBox(b) => b.add_as_element(&mut objects),
        }
    }

    // the BHVNode creator will retuned a Box<Hittable>, either are BHVNode or a Element
    log::info!("Creating the BHV node tree");
    let end = objects.len();
    let bhv_tree: Box<dyn Hittable + Sync> = BHVNode::new(&mut objects, 0, end);

    // create a 1-d vector holding all the pixels, and split into bands for parallel rendering
    let mut pixels =
        vec![Color::new(0.0, 0.0, 0.0); (config.img_width * config.img_height) as usize];
    let bands: Vec<(usize, &mut [Color])> = pixels
        .chunks_mut(config.img_width as usize)
        .enumerate()
        .collect();

    // start the actual render
    log::info!("Starting render!");
    let pb: ProgressBar = ProgressBar::new(config.img_height as u64);

    // use Rayon parallel iterator to iterate over the bands and render per line
    let start: Instant = Instant::now();
    bands.into_par_iter().for_each(|(i, band)| {
        render_line(i as i64, band, &camera, &bhv_tree, &config);
        pb.inc(1);
    });
    
    log::info!("Render finished in {}ms", start.elapsed().as_millis());
    
    // return the pixels
    Ok(pixels)

} // fn render

// fn render_line
// takes a band of pixels (a single horizontal line along a fixed y coordinate) and renders the pixels on that band
pub fn render_line(
    y: i64,
    band: &mut [Color],
    camera: &Camera,
    bhv_tree: &Box<dyn elements::Hittable + Sync>,
    config: &Config,
) {
    let mut rng = SmallRng::seed_from_u64(y as u64);

    for x in 0..config.img_width as usize {
        // start with black color
        let mut color: Color = Color::new(0.0, 0.0, 0.0);

        // loop through all the anti aliasing samples
        for _i in 0..config.samples {
            // get multiple rays for anti alliasing, and add the colors
            let ray: Ray = camera.get_prime_ray(x as i64, y, &mut rng);
            color += ray_color(&bhv_tree, &config, &ray, config.max_depth, &mut rng);
        }

        // set pixel color, but first divide by the number of samples to get the average and return
        band[x] = color.divide_by_samples(config.samples).clamp().linear_to_gamma(2.2);
    }
} // fn render_line

// fn ray_color
// finds the color of the ray, by checking for hits and using the material of that hit + scattered/reflected rays to establish the color
// limit the number of rays it will scatter/reflect by setting depth (e.g. to 32)
fn ray_color(
    bhv_tree: &Box<dyn elements::Hittable + Sync>,
    config: &Config,
    ray: &Ray,
    depth: usize,
    rng: &mut impl Rng,
) -> Color {
    if depth == 0 {
        // we ran out of depth iterations, so we return black
        return Color::new(0.0, 0.0, 0.0);
    }

    // we trace the ray to see what it will hit. uses the BHV_tree to figure our what the first hit is
    let trace = bhv_tree.hit(&ray, &mut Interval::new(0.001, f64::MAX));
    //let trace = scene.trace(&ray);

    match trace {
        Some(hit) => {
            let mut color_from_scatter = Color::new(0.0, 0.0, 0.0);
            // we hit something
            // we start with scattered rays (we assume every object has scattered rays, although in some materials (like metal) its actually a reflected or refracted (glass) ray)
            match hit.material.scatter(ray, &hit, rng) {
                Some((scattered_ray, pdf_option, albedo)) => {
                    // see if there is a ray returned
                    match scattered_ray {
                        Some(sr) => {
                            // there is a scattered ray, so lets get the color of that ray
                            // call the ray_color function again, now with the reflected ray but decrease the depth by 1 so that we dont run into an infinite loop
                            let target_color: Color =
                                ray_color(&bhv_tree, &config, &sr, depth - 1, rng);

                            // get the 
                            let scattering_pdf: f64 = hit.material.scattering_pdf(ray, &hit, &sr);
                            
                            // if there was a pdf included in the scattered ray, apply that one, else use the scattering pdf as the basis (cancelling each other out)
                            let pdf: f64 = if let Some(p) = pdf_option {
                                p
                            } else {
                                scattering_pdf
                            };

                            // return the color, by applying the albedo to the color of the scattered ray (albedo is here defined the amount of color not absorbed)
                            color_from_scatter = Color::new(
                                (albedo.r * target_color.r * scattering_pdf) / pdf,
                                (albedo.g * target_color.g * scattering_pdf) / pdf,
                                (albedo.b * target_color.b * scattering_pdf) / pdf,
                            );
                        }
                        None => {
                            // there is no ray, so just return the albedo of the hittable object (this could be a light?)
                            // albedo
                        }
                    }
                }
                None => {
                    // no scattered ray
                    // return black
                    //return Color::new(0.0, 0.0, 0.0);
                }
            }

            //now its time for emission
            let mut color_from_emmission = Color::new(0.0, 0.0, 0.0);

            match hit.material.emitted(ray, &hit) {
                Some(c) => {
                    color_from_emmission = c;
                }
                _ => {
                    //
                }
            }

            // return the emmission and the scatter color
            color_from_scatter + color_from_emmission
        }
        None => {
            // we did not hit anything, so we return the color of the sky but with a little gradient
            //let a = 0.5 * (ray.direction.y() + 1.0);
            //return Color::new(0.9, 0.9, 1.0) * (1.0 - a) + config.sky_color * a;
            Color::new(0.0, 0.0, 0.0)
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
        //render(scene, config)?;
        Ok(())
    }
}
