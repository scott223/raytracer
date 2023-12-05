use raytracer::render::RenderIntegrator;

// App main function
// load the JSON files for the Scene and for the Configuration, and calls the main raytracer function render from lib.rs
fn main() {
    dotenv::dotenv().ok();
    env_logger::init();

    log::info!(
        "Program started: loading scene and config from JSON files located in /input and render will be located in /render. Current directory is {:?}",
        std::env::current_dir().unwrap()
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

#[cfg(test)]
mod tests {
    use raytracer::{
        bvh::BVHSplitMethod,
        elements::{JSONElement, JSONObj, Rotate, Scale, Transpose},
        linalg::Vec3,
        materials::Lambertian,
        render::{Color, Config, JSONCamera, JSONScene, RenderIntegrator},
    };

    #[test_log::test]
    fn render() {
        let config: Config = Config {
            img_width: 600.,
            img_height: 600.,
            sample_batch_size: 12,
            max_sample_batches: 1,
            min_sample_batches: 1,
            max_depth: 8,
            sky_color: Color::new(0.5, 0.5, 0.5),
            pixel_radius: 2.0,
            bvh_split_method: Some(BVHSplitMethod::Mid),
            gamma_correction: 2.2,
        };
        let json_camera: JSONCamera = JSONCamera {
            camera_center: Vec3::new(278., 278., -800.),
            camera_look_at: Vec3::new(278., 278., 0.),
            camera_fov_vertical: 40.,
            camera_defocus_angle: 0.,
            camera_focus_dist: 800.,
        };

        let dragon = JSONElement::JSONObj(JSONObj {
            filepath: "input/obj/dragon.obj".to_string(),
            transpose: Some(Transpose {
                x: 220.,
                y: 0.,
                z: 200.,
            }),
            rotate: Some(Rotate {
                theta_x: 0.,
                theta_y: -3.2,
                theta_z: 0.,
            }),
            scale: Some(Scale {
                x: 20.,
                y: 20.,
                z: 20.,
            }),
            material: raytracer::materials::Material::Lambertian(Lambertian::new(Color::new(
                0.7, 0.7, 0.7,
            ))),
        });

        let mut elements: Vec<JSONElement> = Vec::new();
        elements.push(dragon);

        let json_scene: JSONScene = JSONScene {
            camera: json_camera,
            elements,
        };

        let mut r = RenderIntegrator::new(json_scene, config);

        // try to execute the main render
        match r.render() {
            Ok(_) => assert!(true),
            Err(_) => assert!(false, "error in render"),
        }
    }
}
