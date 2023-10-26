use pdf::{HittablePDF, MixedPDF};
// Rust imports
use rand::rngs::SmallRng;
use rand::{SeedableRng};

use indicatif::ProgressBar;

use rayon::prelude::*;
use std::error::Error;

use std::time::Instant;

// Raytracer imports
pub mod config;

mod aabb;
mod bhv;
mod camera;
pub mod color;
pub mod elements;
pub mod interval;
mod mat4;
pub mod materials;
mod onb;
pub mod ray;
pub mod vec3;
mod pdf;

use bhv::BHVNode;
use camera::Camera;
use color::Color;
use config::{Config, JSONScene};
use elements::{Element, Hittable};
use interval::Interval;
use materials::{Emmits, Scatterable};
use ray::Ray;

use crate::elements::JSONElement;

use crate::pdf::{PDFTrait, PDF, CosinePDF};

// fn render
// the main render function that sets up the camera, creates an 1d vector for the pixels, splits it into bands, calls the band render function and writes to an image file
// applies parallel rendering using Rayon
pub fn render(json_scene: JSONScene, config: Config) -> Result<Vec<Color>, Box<dyn Error>> {
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

    let mut lights: Vec<Element> = Vec::new();

    for json_element in json_scene.lights {
        match json_element {
            JSONElement::JSONQuad(q) => q.add_as_element(&mut lights),
            JSONElement::JSONSphere(s) => s.add_as_element(&mut lights),
            JSONElement::JSONTriangle(t) => t.add_as_element(&mut lights),
            JSONElement::JSONBox(b) => b.add_as_element(&mut lights),
        }

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
        render_line(i as i64, band, &camera, &bhv_tree, &lights, &config);
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
    lights: &Vec<Element>,
    config: &Config,
) {
    let mut rng: SmallRng = SmallRng::seed_from_u64(y as u64);

    for x in 0..config.img_width as usize {
        // start with black color
        let mut color: Color = Color::new(0.0, 0.0, 0.0);

        // sets a pixel to follow and print detailed logs
        let follow_coords: [usize; 2] = [75, 200];

        // loop through all the anti aliasing samples
        for i in 0..config.samples {
            // get multiple rays for anti alliasing, and add the colors
            let ray: Ray = camera.get_prime_ray(x as i64, y, &mut rng);
            let follow = if x == follow_coords[0] && y as usize == follow_coords[1] { true } else { false };
            if follow  { log::info!("sample: {}", i); }

            color += ray_color(&bhv_tree, &lights, &config, &ray, config.max_depth, &mut rng, follow);
        }

        if color.has_nan() {
            log::error!("color has NaN");
        }

        // set pixel color, but first divide by the number of samples to get the average and return
        band[x] = color
            .divide_by_samples(config.samples)
            .clamp()
            .linear_to_gamma(2.2);
    }
} // fn render_line

// fn ray_color
// finds the color of the ray, by checking for hits and using the material of that hit + scattered/reflected rays to establish the color
// limit the number of rays it will scatter/reflect by setting depth (e.g. to 32)
fn ray_color(
    bhv_tree: &Box<dyn elements::Hittable + Sync>,
    lights: &Vec<Element>,
    config: &Config,
    ray: &Ray,
    depth: usize,
    rng: &mut SmallRng,
    follow: bool,
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
            // we hit something

            if follow {
                log::info!("follow pixel: we hit something");
            }

            // we hit a light, so we return the color from emmission and end here
            if let Some(c) = hit.material.emitted(ray, &hit) {
                
                if follow {
                    log::info!("emmitting, returning color c {:?}", c);
                }
                
                return c;

            }

            // we see if we get a scatter from the material
            if let Some((Some(sscattered), Some(pdf_val), attenuation)) = hit.material.scatter(ray, &hit, rng) {

                //  scattered rays (we assume every object has scattered rays, although in some materials (like metal) its actually a reflected or refracted (glass) ray)

                let light_pdf: PDF = PDF::HittablePDF(HittablePDF::new(hit.point, lights[0]));
                let cosine_pdf: PDF = PDF::CosinePDF(CosinePDF::new(hit.normal));
                let mixed_pdf: PDF = PDF::MixedPDF(MixedPDF::new(hit.point, light_pdf, cosine_pdf));

                // generate the direction for the new scattered ray based on the PDF
                let scattered = Ray::new(hit.point, mixed_pdf.generate(rng));

                // get the pdf value
                let pdf_val = mixed_pdf.value(scattered.direction);


                // there is a scattered ray, so lets get the color of that ray
                // call the ray_color function again, now with the reflected ray but decrease the depth by 1 so that we dont run into an infinite loop
                let target_color: Color =
                    ray_color(&bhv_tree, &lights, &config, &scattered, depth - 1, rng, follow);

                // get the pdf for that spot on that material
                let scattering_pdf: f64 =
                    hit.material.scattering_pdf(ray, &hit, &scattered);

                // return the color, by applying the albedo to the color of the scattered ray (albedo is here defined the amount of color not absorbed)
                // also apply the sample importance weighting
                let color_from_scatter = (attenuation * target_color * scattering_pdf)/ pdf_val;

                // if we track this pixel, print out extra info
                if follow {
                    log::info!("color_from_scatter: {:?}, target color: {:?}, attentuation: {:?}, scattering_pdf: {:?}, pdf_val: {:?}", color_from_scatter, target_color, attenuation, scattering_pdf, pdf_val);
                 }

                 if color_from_scatter.has_nan() { log::error!("color has Nan!"); }

                return color_from_scatter;


            }

            // no emmission, nor scatter. just return black
            Color::new(0.0, 0.0, 0.0)

        },
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
