use crate::{vec3::Vec3, ray::Ray};
use rand::{Rng};

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
impl Camera {
    pub fn new(img_width: f64, sensor_ratio: f64) -> Self {
        let camera_center = Vec3::new(0.0,0.0,0.0);
        let focal_length: f64 = 5.0;
        let vv: f64 = 2.0; // viewport height 
        let vu: f64 = vv * sensor_ratio; // calculating viewport width

        // Calculate the vectors across the horizontal and down the vertical viewport edges. this is basically the coordinate of the two corners
        let viewport_u: Vec3 = Vec3::new(vu, 0.0, 0.0);
        let viewport_v: Vec3 = Vec3::new(0.0, -vv, 0.0);

        // spacing between each pixel
        let pixel_delta_u: Vec3 = viewport_u / (img_width as f64);
        let pixel_delta_v: Vec3 = viewport_v / ((img_width/sensor_ratio) as f64);   

        // now we position the viewport in front of the camera, taking into account the camera position, the focal length, and the viewport size (we devide those by 2 to put the camera in the middle)
        let viewport_upper_left = camera_center - Vec3::new(0.0,0.0,focal_length) - viewport_u/2.0 - viewport_v/2.0;
        let pixel00_loc = viewport_upper_left + (pixel_delta_u + pixel_delta_v) * 0.5;     

        log::info!("sensor ratio: {:.3}, vv: {:.3}, vu {:.3}, pixel_delta_u {}, pixel_delta_v {}", sensor_ratio, vv, vu, pixel_delta_u, pixel_delta_v);

        Camera{
            viewport_u: viewport_u,
            viewport_v: viewport_v,
            focal_length: focal_length,
            camera_center: camera_center,
            viewport_upper_left: viewport_upper_left,
            pixel_delta_u: pixel_delta_u,
            pixel_delta_v: pixel_delta_v,
            pixel00_loc: pixel00_loc,
        }

    }

    // get a prime ray from the camera, add a stochastic x and y component for anti allisasing
    pub fn get_prime_ray(self, x: i64, y: i64) -> Ray {
            let pixel_loc = self.pixel00_loc + (self.pixel_delta_u * x as f64) + (self.pixel_delta_v * y as f64);
            let pixel_sample = pixel_loc + self.pixel_sample_square();
            let ray_direction = pixel_sample - self.camera_center;
            
            let ray = Ray::new(self.camera_center, ray_direction);
            // log::trace!("Ray created for pixel ({}.{}): {}", x, y, ray);

            return ray;
    }

    // generate a random point inside the box of -0.5 and +0.5 units * pixel size around the pixel center
    pub fn pixel_sample_square(self) -> Vec3 {
        let mut rng = rand::thread_rng();

        let n1: f64 = rng.gen_range(-0.5..0.5);
        let n2: f64 = rng.gen_range(-0.5..0.5);

        let point = Vec3::new(n1*self.pixel_delta_u.x(), n2*self.pixel_delta_v.y(), 0.0);

        return point;
    }
}

#[cfg(test)]
mod tests {
    use crate::vec3::Vec3;
    use crate::camera::Camera;
    use assert_approx_eq::assert_approx_eq;

    #[test_log::test]
    fn test_create_camera() {
        // TODO lots of camera tests
        let camera: Camera = Camera::new(800.0, 9.0/16.0);

        assert_approx_eq!(camera.camera_center, Vec3::new(0.0, 0.0, 0.0));
    }
}