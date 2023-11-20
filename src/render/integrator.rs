use std::{
    error::Error,
    fs::File,
    io::{BufReader, BufWriter},
    path::Path,
    time::Instant,
};

use indicatif::ProgressBar;
use rand::{rngs::SmallRng, SeedableRng};
use rayon::prelude::*;

use crate::{
    bvh::BVH_SAH,
    elements::{Element, JSONElement},
    materials::{Emmits, Reflects, Refracts, Scatterable},
    render::camera::Camera,
    render::pdf::{HittablePDF, MixedPDF, PDFTrait, Pdf},
    render::ray::Ray,
    render::Color,
    render::Interval,
    render::{Config, JSONScene},
};

#[derive(Debug, Clone)]
pub struct RenderIntegrator {
    json_scene: JSONScene,
    config: Config,
    pixels: Vec<Color>,
}

impl RenderIntegrator {
    pub fn new(json_scene: JSONScene, config: Config) -> Self {
        // create a 1-d vector holding all the pixels, and
        let pixels =
            vec![Color::new(0.0, 0.0, 0.0); (config.img_width * config.img_height) as usize];

        RenderIntegrator {
            json_scene,
            config,
            pixels,
        }
    }

    pub fn new_from_json(
        scene_path: &str,
        config_path: &str,
    ) -> Result<RenderIntegrator, Box<dyn Error>> {
        // Open the Scene file
        let scene_file: File = File::open(scene_path)?;
        let scene_reader: BufReader<File> = BufReader::new(scene_file);

        // Read the JSON contents of the file as an instance of `Scene`.
        let scene: JSONScene = serde_json::from_reader(scene_reader)?;

        // Open the Config file
        let config_file: File = File::open(config_path)?;
        let config_reader: BufReader<File> = BufReader::new(config_file);

        // Read the JSON contents of the file as an instance of `Config`.
        let config: Config = serde_json::from_reader(config_reader)?;

        return Ok(RenderIntegrator::new(scene, config));
    }

    pub fn save_to_png(&self, png_path: &str) -> Result<(), Box<dyn Error>> {
        let mut raw_pixels: Vec<u8> = Vec::new();

        for p in &self.pixels {
            let rgb = p.to_rgb();
            raw_pixels.push(rgb.r);
            raw_pixels.push(rgb.g);
            raw_pixels.push(rgb.b);
        }

        let path = Path::new(png_path);
        let file = File::create(path).unwrap();
        let w = &mut BufWriter::new(file);

        let mut encoder = png::Encoder::new(
            w,
            self.config.img_width as u32,
            self.config.img_height as u32,
        ); // Width x heigth
        encoder.set_color(png::ColorType::Rgb);
        encoder.set_depth(png::BitDepth::Eight);
        encoder.set_source_gamma(png::ScaledFloat::new(1.0 / 2.2)); // 1.0 / 2.2, unscaled, but rounded
        let source_chromaticities = png::SourceChromaticities::new(
            // Using unscaled instantiation here
            (0.31270, 0.32900),
            (0.64000, 0.33000),
            (0.30000, 0.60000),
            (0.15000, 0.06000),
        );
        encoder.set_source_chromaticities(source_chromaticities);
        let mut writer = encoder.write_header().unwrap();

        let data = raw_pixels; // An array containing a RGB sequence
        writer.write_image_data(&data).unwrap(); // Save

        Ok(())
    }

    // fn render
    // the main render function that sets up the camera, creates an 1d vector for the pixels, splits it into bands, calls the band render function and writes to an image file
    // applies parallel rendering using Rayon
    pub fn render(&mut self) -> Result<(), Box<dyn Error>> {
        // create a new camera object
        let camera: Camera = Camera::new(&self.config, &self.json_scene.camera);

        // we loop through the JSON Scene, and create a new actual object for each simplified JSON object
        // this allows us to pre-calculate various things like the bbox or commonly used numbers in the cconstructor
        let mut objects: Vec<Element> = Vec::new();

        for json_element in &self.json_scene.elements {
            match json_element {
                JSONElement::JSONQuad(q) => q.add_as_element(&mut objects),
                JSONElement::JSONSphere(s) => s.add_as_element(&mut objects),
                JSONElement::JSONTriangle(t) => t.add_as_element(&mut objects),
                JSONElement::JSONBox(b) => b.add_as_element(&mut objects),
                JSONElement::JSONObj(o) => o.add_as_element(&mut objects),
            }
        }

        // the BHVNode creator will return a Box<Hittable>, either are BHVNode (if there are more than 1 objects) or an Element
        log::info!("Creating the BVH node tree, with {} objects", objects.len());

        let bvh_tree_sah = BVH_SAH::build(&objects);

        // create a refrence to the elements that are marked as an attractor
        let attractors: Vec<&Element> = objects.iter().filter(|e| e.is_attractor()).collect();

        // split into bands for parallel rendering

        let bands: Vec<(usize, &mut [Color])> = self
            .pixels
            .chunks_mut(self.config.img_width as usize)
            .enumerate()
            .collect();

        // start the actual render
        log::info!("Starting render");
        let start: Instant = Instant::now();
        let pb: ProgressBar = ProgressBar::new(self.config.img_height as u64);

        // use Rayon parallel iterator to iterate over the bands and render per line
        bands.into_par_iter().for_each(|(i, band)| {
            Self::render_line(
                i as i64,
                band,
                &camera,
                &bvh_tree_sah,
                &attractors,
                &self.config,
            );
            pb.inc(1);
        });

        log::info!("Render finished in {}ms", start.elapsed().as_millis());

        // return an OK
        Ok(())
    }

    // fn render_line
    // takes a band of pixels (a single horizontal line along a fixed y coordinate) and renders the pixels on that band
    fn render_line(
        y: i64,
        band: &mut [Color],
        camera: &Camera,
        bvh_tree_sah: &BVH_SAH,
        lights: &Vec<&Element>,
        config: &Config,
    ) {
        let mut rng: SmallRng = SmallRng::seed_from_u64(y as u64);

        for (x, band_item) in band.iter_mut().enumerate().take(config.img_width as usize) {
            // start with black color
            let mut color: Color = Color::new(0.0, 0.0, 0.0);

            // sets a pixel to follow and print detailed logs
            let follow_coords: [usize; 2] = [40, 120];

            // loop through all the anti aliasing samples
            for i in 0..config.samples {
                // get multiple rays for anti alliasing, and add the colors
                let ray: Ray = camera.get_prime_ray(x as i64, y, &mut rng);
                let follow: bool = x == follow_coords[0] && y as usize == follow_coords[1];
                //let follow = false;

                if follow {
                    log::info!("sample: {}", i);
                }

                color += Self::ray_color(
                    bvh_tree_sah,
                    lights,
                    &ray,
                    config.max_depth,
                    &mut rng,
                    follow,
                );
            }

            // set pixel color, but first divide by the number of samples to get the average and return
            *band_item = color
                .divide_by_samples(config.samples)
                .clamp()
                .linear_to_gamma(2.2);
        }
    }

    // fn ray_color
    // finds the color of the ray, by checking for hits and using the material of that hit + scattered/reflected rays to establish the color
    // limit the number of rays it will scatter/reflect by setting depth (e.g. to 32)
    fn ray_color(
        bvh_tree_sah: &BVH_SAH,
        attractors: &Vec<&Element>,
        ray: &Ray,
        depth: usize,
        rng: &mut SmallRng,
        follow: bool,
    ) -> Color {
        if depth == 0 {
            log::trace!("Ran out of depth");
            // we ran out of depth iterations, so we return black
            return Color::new(0.0, 0.0, 0.0);
        }

        // we trace the ray to see what it will hit. uses the BHV_tree to figure our what the first hit is
        let trace = bvh_tree_sah.hit(ray, &mut Interval::new(0.0001, f64::MAX), follow);

        match trace {
            Some(hit) => {
                // we hit something
                if follow {
                    log::info!("follow pixel: we hit something");
                }

                // we hit a light, so we return the color from emmission and end here
                if let Some(color_from_emission) = hit.material.emitted(ray, &hit) {
                    if follow {
                        log::info!("emmitting, returning color c {:?}", color_from_emission);
                    }

                    return color_from_emission;
                }

                // we see if we get a scatter from the material
                if let Some(scatter) = hit.material.scatter(ray, &hit, rng) {
                    // generate the direction for the new scattered ray based on the PDF
                    // check if there are attractors, else just use the scatter pdf
                    let (scattered_ray, pdf_val) = if attractors.len() > 0 {
                        let hittable_pdf =
                            Pdf::HittablePDF(HittablePDF::new(hit.point, attractors));
                        let mixed_pdf = Pdf::MixedPDF(MixedPDF::new(
                            hit.point,
                            0.5,
                            &hittable_pdf,
                            &scatter.pdf,
                        ));
                        let scattered_ray = Ray::new(hit.point, mixed_pdf.generate(rng));
                        (scattered_ray, mixed_pdf.value(scattered_ray.direction))
                    } else {
                        let scattered_ray = Ray::new(hit.point, scatter.pdf.generate(rng));
                        (scattered_ray, scatter.pdf.value(scattered_ray.direction))
                    };

                    // there is a scattered ray, so lets get the color of that ray
                    // call the ray_color function again, now with the reflected ray but decrease the depth by 1 so that we dont run into an infinite loop
                    let target_color: Color = Self::ray_color(
                        bvh_tree_sah,
                        attractors,
                        &scattered_ray,
                        depth - 1,
                        rng,
                        follow,
                    );

                    // get the pdf for that spot on that material
                    let scattering_pdf: f64 = scatter.pdf.value(scattered_ray.direction);

                    // return the color, by applying the albedo to the color of the scattered ray (albedo is here defined the amount of color not absorbed)
                    // also apply the sample importance weighting
                    let mut color_from_scatter =
                        (scatter.attenuation * target_color * scattering_pdf) / pdf_val;

                    // there could be some NaNs as pdf_val = 0 is a div by zero. took the lazy route and just filtere them out
                    if color_from_scatter.has_nan() {
                        color_from_scatter = Color::new(0., 0., 0.);
                    }

                    // if we track this pixel, print out extra info
                    if follow {
                        log::info!("color_from_scatter: {}, target color: {}, attentuation: {}, scattering_pdf: {:.3}, pdf_val: {:.3}", color_from_scatter, target_color, scatter.attenuation, scattering_pdf, pdf_val);
                    }

                    return color_from_scatter;
                }

                // did we get a reflection?
                if let Some(reflect) = hit.material.reflect(ray, &hit, rng) {
                    if let Some(reflected_ray) = reflect.ray {
                        // there is refleced ray, so lets follow that one
                        let target_color: Color = Self::ray_color(
                            bvh_tree_sah,
                            attractors,
                            &reflected_ray,
                            depth - 1,
                            rng,
                            follow,
                        );
                        return reflect.attenuation * target_color;
                    } else {
                        //no reflected ray, just return the attentuation (this is now set to 0, 0, 0)
                        return reflect.attenuation;
                    }
                }

                // and finally refraction
                if let Some(refract) = hit.material.refract(ray, &hit, rng) {
                    let target_color: Color = Self::ray_color(
                        bvh_tree_sah,
                        attractors,
                        &refract.ray,
                        depth - 1,
                        rng,
                        follow,
                    );

                    return refract.attenuation * target_color;
                }

                // no emmission, nor scatter. just return black
                Color::new(0.0, 0.0, 0.0)
            }
            None => {
                // we did not hit anything, so we return black for now
                Color::new(0.0, 0.0, 0.0)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    //use crate::config::Config;
    //use crate::config::Scene;

    use std::error::Error;

    #[test_log::test]
    fn test_render_full_scene() -> Result<(), Box<dyn Error>> {
        //TODO: test and default for config and scene

        //let _config: Config = Config::default();
        //let _scene: Scene = Scene::default();
        //render(scene, config)?;
        Ok(())
    }
}
