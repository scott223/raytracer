use crate::{config::Config, ray::Ray, vec3::Vec3};
use rand::Rng;

#[derive(Debug, Clone, Copy)]
pub struct Camera {
    pub camera_center: Vec3,
    pub focal_length: f64,
    pub viewport_u: Vec3,
    pub viewport_v: Vec3,
    pub viewport_upper_left: Vec3,
    pub pixel_delta_u: Vec3,
    pub pixel_delta_v: Vec3,
    pub pixel00_loc: Vec3,
}

// Create a new Camera with a origin and a direction (normalized)
// this function setts up the viewport, and calculates the location of the upper left pixel (00)
impl Camera {
    pub fn new(config: &Config) -> Self {
        let camera_center: Vec3 = config.camera_center;
        let focal_length: f64 = config.focal_length;

        let ratio = config.img_width / config.img_height;

        let vv: f64 = 2.0; // viewport height
        let vu: f64 = vv * ratio; // calculating viewport width

        // Calculate the vectors across the horizontal and down the vertical viewport edges. this is basically the coordinate of the two corners
        let viewport_u: Vec3 = Vec3::new(vu, 0.0, 0.0);
        let viewport_v: Vec3 = Vec3::new(0.0, -vv, 0.0);

        // spacing between each pixel
        let pixel_delta_u: Vec3 = viewport_u / (config.img_width as f64);
        let pixel_delta_v: Vec3 = viewport_v / ((config.img_width / ratio) as f64);

        // now we position the viewport in front of the camera, taking into account the camera position, the focal length, and the viewport size (we devide those by 2 to put the camera in the middle)
        let viewport_upper_left =
            camera_center - Vec3::new(0.0, 0.0, focal_length) - viewport_u / 2.0 - viewport_v / 2.0;
        let pixel00_loc = viewport_upper_left + (pixel_delta_u + pixel_delta_v) * 0.5;

        log::info!(
            "sensor ratio: {:.3}, vv: {:.3}, vu {:.3}, pixel_delta_u {}, pixel_delta_v {}",
            ratio,
            vv,
            vu,
            pixel_delta_u,
            pixel_delta_v
        );

        Camera {
            viewport_u,
            viewport_v,
            focal_length,
            camera_center,
            viewport_upper_left,
            pixel_delta_u,
            pixel_delta_v,
            pixel00_loc,
        }
    }

    // get a prime ray from the camera, add a stochastic x and y component for anti allisasing
    pub fn get_prime_ray(self, x: i64, y: i64) -> Ray {
        let pixel_loc =
            self.pixel00_loc + (self.pixel_delta_u * x as f64) + (self.pixel_delta_v * y as f64);
        let pixel_sample = pixel_loc + self.pixel_sample_square();
        let ray_direction = pixel_sample - self.camera_center;

        Ray::new(self.camera_center, ray_direction)
    }

    // generate a random point inside the box of -0.5 and +0.5 units * pixel size around the pixel center
    // TODO: this function now seems to take a lot of the time in a thread, so need to explore how to do this faster
    pub fn pixel_sample_square(self) -> Vec3 {
        let mut rng = rand::thread_rng();

        let n1: f64 = rng.gen_range(-0.5..0.5);
        let n2: f64 = rng.gen_range(-0.5..0.5);

        let point = Vec3::new(
            n1 * self.pixel_delta_u.x(),
            n2 * self.pixel_delta_v.y(),
            0.0,
        );

        return point;
    }
}

#[cfg(test)]
mod tests {
    use crate::camera::Camera;
    use crate::config::Config;
    use crate::vec3::Vec3;
    use assert_approx_eq::assert_approx_eq;

    #[test_log::test]
    fn test_create_camera() {
        let config: Config = Config {
            img_width: 1024.0,
            img_height: 576.0,
            camera_center: Vec3::new(0.0, 0.0, 0.0),
            focal_length: 5.0,
            samples: 1,
            max_depth: 32,
        };

        let camera: Camera = Camera::new(&config);

        assert_approx_eq!(camera.camera_center, Vec3::new(0.0, 0.0, 0.0));
        assert_approx_eq!(camera.focal_length, 5.0);

        assert_approx_eq!(camera.viewport_v, Vec3::new(0.0, -2.0, 0.0));
        assert_approx_eq!(camera.viewport_u, Vec3::new(3.555555555555555, 0.0, 0.0));
        assert_approx_eq!(camera.viewport_upper_left, Vec3::new(-1.77777777777, 1.0, -5.0));
        assert_approx_eq!(camera.pixel_delta_u, Vec3::new(0.00347222, 0.0, 0.0));
        assert_approx_eq!(camera.pixel_delta_v, Vec3::new(0.0, -0.0034722222, 0.0));
        assert_approx_eq!(camera.pixel00_loc, Vec3::new(-1.77604166666666, 0.998263888, -5.0));
    }
}