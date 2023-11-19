use dotenv;
use std::env;

use raytracer::render::RenderIntegrator;

// App main function
// load the JSON files for the Scene and for the Configuration, and calls the main raytracer function render from lib.rs
fn main() {
    dotenv::dotenv().ok();
    env_logger::init();

    log::info!(
        "Program started: loading scene and config from JSON files located in /input and render will be located in /render. Current directory is {:?}",
        env::current_dir().unwrap()
    );

    //load the JSON into the integrator
    match RenderIntegrator::new_from_json("input/scene.json", "input/config.json") {
        Ok(mut r) => {
            // succes, so execute the main render
            match r.render() {
                Ok(()) => {
                    // success, so save to png
                    match r.save_to_png("renders/render.png") {
                        Ok(()) => {
                            log::info!("File saved, closing program");
                        }
                        Err(e) => {
                            log::error!("Error in saving file: {}", e);
                        }
                    }
                }
                Err(e) => {
                    log::error!("Error in render: {}", e);
                }
            };
        }
        Err(e) => {
            log::error!("Error in reading JSON: {}", e);
        }
    }
} //fn main
