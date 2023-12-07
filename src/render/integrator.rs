use std::{
    error::Error,
    fmt::Write,
    fs::File,
    io::{BufReader, BufWriter},
    path::Path,
};

use indicatif::{ProgressBar, ProgressState, ProgressStyle};
use rand::{rngs::SmallRng, SeedableRng};
use rayon::prelude::*;

use crate::{
    bvh::{BVHSplitMethod, BVH_SAH},
    elements::{Element, JSONElement},
    materials::{Emmits, Reflects, Refracts, Scatterable},
    render::camera::Camera,
    render::pdf::{HittablePDF, MixedPDF, PDFTrait, Pdf},
    render::ray::Ray,
    render::Color,
    render::Interval,
    render::{Config, JSONScene},
};

use sobol_burley::sample;

use super::{filter::{Filter, MitchellNetravali}, Stats};

#[derive(Debug, Clone)]
pub struct RenderIntegrator {
    json_scene: JSONScene,
    config: Config,
    pixels: Vec<Color>,
    sample_pixels: Vec<Color>,
    filter: Filter,
}

impl RenderIntegrator {
    pub fn new(json_scene: JSONScene, config: Config) -> Self {
        // create a 1-d vector holding all the pixels, and
        let pixels =
            vec![Color::new(0.0, 0.0, 0.0); (config.img_width * config.img_height) as usize];

        let sample_pixels =
            vec![Color::new(0.0, 0.0, 0.0); (config.img_width * config.img_height) as usize];

        let filter: Filter = Filter::MitchellNetravali(MitchellNetravali::new(
            config.pixel_radius,
            config.pixel_radius,
            1. / 3.,
            1. / 3.,
        ));

        RenderIntegrator {
            json_scene,
            config,
            pixels,
            sample_pixels,
            filter,
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

        Ok(RenderIntegrator::new(scene, config))
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

        // Sample image

        let mut raw_pixels: Vec<u8> = Vec::new();

        for p in &self.sample_pixels {
            let rgb = p.to_rgb();
            raw_pixels.push(rgb.r);
            raw_pixels.push(rgb.g);
            raw_pixels.push(rgb.b);
        }

        let path = Path::new("samples.png");
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

        // see if we have a split method in the config, else pick a default
        let split_method = if self.config.bvh_split_method.is_some() {
            self.config.bvh_split_method.unwrap()
        } else {
            BVHSplitMethod::Mid
        };

        let bvh_tree_sah = BVH_SAH::build(&objects, split_method);

        // create a refrence to the elements that are marked as an attractor
        let attractors: Vec<&Element> = objects.iter().filter(|e| e.is_attractor()).collect();

        // split into bands for parallel rendering

        let bands: Vec<(usize, &mut [Color])> = self
            .pixels
            .chunks_mut(self.config.img_width as usize)
            .enumerate()
            .collect();

        let sample_bands: Vec<(usize, &mut [Color])> = self
            .sample_pixels
            .chunks_mut(self.config.img_width as usize)
            .enumerate()
            .collect();

        // start the actual render
        log::info!("Starting render");
        let pb: ProgressBar = ProgressBar::new(self.config.img_height as u64);

        pb.set_style(
            ProgressStyle::with_template(
                "{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({eta})",
            )
            .unwrap()
            .with_key("eta", |state: &ProgressState, w: &mut dyn Write| {
                write!(w, "{:.1}s", state.eta().as_secs_f64()).unwrap()
            })
            .progress_chars("#>-"),
        );

        //let mut sample_pixels: Vec<Color> = Vec::with_capacity((self.config.img_height * self.config.img_width) as usize);

        // use Rayon parallel iterator to iterate over the bands and render per line
        bands
            .into_par_iter()
            .zip(sample_bands.into_par_iter())
            .enumerate()
            .for_each(|(i, (band, sample_band))| {
                Self::render_line(
                    &self.filter,
                    i,
                    band.1,
                    sample_band.1,
                    &camera,
                    &bvh_tree_sah,
                    &attractors,
                    &self.config,
                );

                pb.inc(1);
            });

        log::info!("Render finished!");

        // return an OK
        Ok(())
    }

    // fn render_line
    // takes a band of pixels (a single horizontal line along a fixed y coordinate) and renders the pixels on that band
    fn render_line(
        filter: &Filter,
        y: usize,
        band: &mut [Color],
        sample_pixels: &mut [Color],
        camera: &Camera,
        bvh_tree_sah: &BVH_SAH,
        lights: &Vec<&Element>,
        config: &Config,
    ) {
        let mut rng: SmallRng = SmallRng::seed_from_u64(y as u64);
        //let mut sample_pixels: Vec<Color> = Vec::with_capacity(config.img_width as usize);

        let max_samples: usize = config.sample_batch_size * config.max_sample_batches;

        for (x, band_item) in band.iter_mut().enumerate() {
            // start with black color

            let mut actual_samples: usize = 0;

            // sets a pixel to follow and print detailed logs
            let _follow_coords: [usize; 2] = [40, 120];
            let pixel_num: usize = x * y;

            let mut stats = Stats::new(max_samples);

            'batch: for b in 0..config.max_sample_batches {
                let mut batch_color: Color = Color::new(0.0, 0.0, 0.0);
                let mut batch_sum_sample_weight: f64 = 0.0;

                // loop through all the anti aliasing samples
                for i in 0..config.sample_batch_size {
                    // use the Sobol sampler to get a sample position relevative to the pixel center
                    let sample_position = vec![
                        ((sample(i as u32, 0, pixel_num as u32) as f64 * config.pixel_radius)
                            - config.pixel_radius / 2.),
                        ((sample(i as u32, 1, pixel_num as u32) as f64 * config.pixel_radius)
                            - config.pixel_radius / 2.),
                    ];

                    let sample_weight = filter.evaluate(sample_position[0], sample_position[1]);

                    // get multiple rays for anti alliasing, and add the colors
                    let ray: Ray = camera.get_prime_ray(x, y, &mut rng, sample_position);

                    //let follow: bool = x == follow_coords[0] && y as usize == follow_coords[1];
                    let follow = false;

                    if follow {
                        log::info!("sample: {}", i);
                    }

                    batch_color += Self::ray_color(
                        bvh_tree_sah,
                        lights,
                        &ray,
                        config.max_depth,
                        &mut rng,
                        follow,
                    ) * sample_weight;

                    batch_sum_sample_weight += sample_weight;
                    actual_samples += 1;
                    
                }

                stats.push(batch_color/batch_sum_sample_weight, batch_sum_sample_weight);

                if b > config.min_sample_batches {
                    //log::info!("var {}, check {}", stats.interval(), 0.20 * stats.mean().illuminance());

                    if stats.interval() < 0.05 * stats.mean().illuminance() {
                        break 'batch;
                    }
                }
            }

            // set pixel color, but first divide by the number of samples to get the average and return
            *band_item = (stats.mean()).clamp().linear_to_gamma(config.gamma_correction);

            let factor: f64 =
                (actual_samples - config.min_sample_batches * config.sample_batch_size) as f64
                    / max_samples as f64;
            //log::info!("sample_factor: {}",factor);
            sample_pixels[x] = Color::new(factor as f64 * 1.0, (1.0 - factor) * 1.0, 0.0);
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
                    let (scattered_ray, pdf_val) = if !attractors.is_empty() {
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
