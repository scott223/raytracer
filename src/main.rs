

fn main() {
    env_logger::init();

    log::info!("Program started");
    let _ = raytracer::render(&raytracer::config::Config::default());
}