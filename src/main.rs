use std::env;
use std::fs::File;
use std::io::BufReader;

// App main function
// load the JSON files for the Scene and for the Configuration, and calls the main raytracer function render from lib.rsß
fn main() {
    env_logger::init();
    log::info!(
        "Program started: loading scene and config from JSON files located in /input. Base directory is {:?}",
        env::current_dir().unwrap()
    );

    // Open the Scene file
    let scene_file =
        File::open("input/scene.json").expect("Error reading input/scene.json, quitting");
    let scene_reader = BufReader::new(scene_file);

    // Read the JSON contents of the file as an instance of `Scene`.
    let scene: raytracer::config::Scene =
        serde_json::from_reader(scene_reader).expect("Error parsing input/scene.json, quitting");

    // Open the Config file
    let config_file =
        File::open("input/config.json").expect("Error reading input/config.json, quitting");
    let config_reader = BufReader::new(config_file);

    // Read the JSON contents of the file as an instance of `Config`.
    let config: raytracer::config::Config =
        serde_json::from_reader(config_reader).expect("Error parsing input/config.json, quitting");

    let _ = raytracer::render(scene, config);
} //fn main
