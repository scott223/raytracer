fn main() {
    env_logger::init();
    log::info!("Program started");
    let _ = raytracer::render();
}