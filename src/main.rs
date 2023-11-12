use std::env;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::Path;

use raytracer::render::integrator::RenderIntegrator;

// App main function
// load the JSON files for the Scene and for the Configuration, and calls the main raytracer function render from lib.rsß
fn main() {
    env_logger::init();

    log::info!(
        "Program started: loading scene and config from JSON files located in /input and render will be located in /render. Base directory is {:?}",
        env::current_dir().unwrap()
    );

    // Open the Scene file
    let scene_file =
        File::open("input/scene.json").expect("Error reading input/scene.json, quitting");
    let scene_reader = BufReader::new(scene_file);

    // Read the JSON contents of the file as an instance of `Scene`.
    let scene: raytracer::config::JSONScene =
        serde_json::from_reader(scene_reader).expect("Error parsing input/scene.json, quitting");

    // Open the Config file
    let config_file =
        File::open("input/config.json").expect("Error reading input/config.json, quitting");
    let config_reader = BufReader::new(config_file);

    // Read the JSON contents of the file as an instance of `Config`.
    let config: raytracer::config::Config =
        serde_json::from_reader(config_reader).expect("Error parsing input/config.json, quitting");

    let r: RenderIntegrator = RenderIntegrator::new(scene, config);

    match r.render() {
        Ok(i) => {
            //i.save("renders/render.png").unwrap();

            let mut raw_pixels: Vec<u8> = Vec::new();

            for p in i {
                let rgb = p.to_rgb();
                raw_pixels.push(rgb.r);
                raw_pixels.push(rgb.g);
                raw_pixels.push(rgb.b);
            }

            let path = Path::new(r"renders/render.png");
            let file = File::create(path).unwrap();
            let w = &mut BufWriter::new(file);

            let mut encoder =
                png::Encoder::new(w, config.img_width as u32, config.img_height as u32); // Width x heigth
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
            log::info!("File saved, closing program");
        }
        Err(e) => {
            panic!("Error in render function: {}", e);
        }
    }
} //fn main
