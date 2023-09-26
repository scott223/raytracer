use rand::Rng;

use crate::{config::Config, ray::Ray, vec3::Vec3};

#[derive(Debug, Clone, Copy)]
pub struct Camera {
    pub camera_center: Vec3, // location of the lense
    pub look_at: Vec3, // focal point that the camera is pointed to
    pub field_of_view: f64, // vertical fov in degrees
    pub relative_up: Vec3, // camera relative "up" position
    pub defocus_angle: f64, // variation angle of rays through each pixel
    pub focus_dist: f64, // distance from camera center to pane of perfect focus
    pub defocus_disk_u: Vec3,  // defocus disk horizontal radius
    pub defocus_disk_v: Vec3,  // defocus disk vertical radius
    pub defocus_radius: f64, // radius of the thin defocus disk
    pub viewport_upper_left: Vec3, // calculated field
    pub pixel_delta_u: Vec3, // calculated field, spacing between each pixel in x
    pub pixel_delta_v: Vec3, // calculated field, spacing between each pixel in y
    pub pixel00_loc: Vec3, // calculated field, location of upper left pixel in world axis

}

// Create a new Camera with a origin and a direction (normalized)
// this function setts up the viewport, and calculates the location of the upper left pixel (00)
impl Camera {
    pub fn new(config: &Config) -> Self {
        let relative_up: Vec3 = Vec3::new(0.0, 1.0, 0.0);

        let camera_center: Vec3 = config.camera_center;
        let look_at = config.camera_look_at;
        let focus_dist = config.camera_focus_dist;

        let field_of_view = config.camera_fov_vertical;

        let ratio = config.img_width / config.img_height;

        let theta = field_of_view.to_radians();
        let h = (theta / 2.0).tan();

        let vv: f64 = 2.0 * h * focus_dist; // viewport height
        let vu: f64 = vv * ratio; //  viewport width

        // Calculate the u,v,w unit basis vectors for the camera coordinate frame.
        let w: Vec3 = (camera_center - look_at).normalized();
        let u: Vec3 = relative_up.cross(&w).normalized();
        let v: Vec3 = w.cross(&u);

        // Calculate the vectors across the horizontal and down the vertical viewport edges. this is basically the coordinate of the two corners
        let viewport_u: Vec3 = u * vu;
        let viewport_v: Vec3 = v * -1.0 * vv;

        // spacing between each pixel
        let pixel_delta_u: Vec3 = viewport_u / (config.img_width as f64);
        let pixel_delta_v: Vec3 = viewport_v / ((config.img_width / ratio) as f64);

        // now we position the viewport in front of the camera, taking into account the camera position, the focal length, and the viewport size (we devide those by 2 to put the camera in the middle)
        let viewport_upper_left: Vec3 =
            camera_center - (w * focus_dist) - viewport_u / 2.0 - viewport_v / 2.0;
        let pixel00_loc: Vec3 = viewport_upper_left + (pixel_delta_u + pixel_delta_v) * 0.5;

        // calculate the defocus disk radius and coordinates
        let defocus_angle = config.camera_defocus_angle;
        let defocus_radius: f64 = focus_dist * ((defocus_angle / 2.0).to_radians().tan());
        let defocus_disk_u: Vec3 = u * defocus_radius;
        let defocus_disk_v: Vec3 = v * defocus_radius;

        let c = Camera {
            relative_up,
            look_at,
            focus_dist,
            camera_center,
            field_of_view,
            defocus_angle,
            defocus_radius,
            defocus_disk_u,
            defocus_disk_v,
            viewport_upper_left,
            pixel_delta_u,
            pixel_delta_v,
            pixel00_loc,
        };

        // TODO beautify display
        log::info!(
            "{:?}", c
        );

        return c; 
    }

    // Get a randomly-sampled camera ray for the pixel at location x,i, originating from
    // the camera defocus disk.
    pub fn get_prime_ray(self, x: i64, y: i64, rng: &mut impl Rng) -> Ray {
        
        // get a random point in the camera defocus disk to be used as an origin, or the camera center if the defocus angle <= 0 (this is for focus blur)
        let p: Vec3 = Vec3::new_random_in_unit_disk(rng);
        let disk_sample: Vec3 = self.camera_center + (self.defocus_disk_u * p.x()) + (self.defocus_disk_v * p.y());
        let ray_origin: Vec3 = if self.defocus_angle <= 0.0 { self.camera_center } else { disk_sample };
        //println!("origin {:?}", ray_origin);

        // establish the the direction of the ray by taking a random sample from the square around the pixel center (this is for anti aliasing)
        let pixel_loc: Vec3 =
            self.pixel00_loc + (self.pixel_delta_u * x as f64) + (self.pixel_delta_v * y as f64);
        let pixel_sample: Vec3 = pixel_loc + self.pixel_sample_square(rng);
        let ray_direction: Vec3 = pixel_sample - ray_origin;

        Ray::new(ray_origin, ray_direction)
    }

    // generate a random point inside the box of -0.5 and +0.5 units * pixel size around the pixel center
    // TODO: this function now seems to take a lot of the time in a thread, so need to explore how to do this faster
    pub fn pixel_sample_square(self, rng: &mut impl Rng) -> Vec3 {

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
    use crate::color::Color;
    use assert_approx_eq::assert_approx_eq;

    #[test_log::test]
    fn test_create_camera() {
        let config: Config = Config {
            img_width: 1024.0,
            img_height: 576.0,
            camera_center: Vec3::new(0.0, 0.0, 0.0),
            camera_look_at: Vec3::new(0.0, 0.0, -5.0),
            camera_fov_vertical: 90.0,
            camera_focus_dist: 5.0,
            camera_defocus_angle: 2.0,
            samples: 1,
            max_depth: 32,
            sky_color: Color::new(3.0 / 255.0, 165.0 / 255.0, 252.0 / 255.0),
        };

        let camera: Camera = Camera::new(&config);

        assert_approx_eq!(camera.camera_center, Vec3::new(0.0, 0.0, 0.0));
        assert_approx_eq!(camera.focus_dist, 5.0);

        // assert_approx_eq!(camera.viewport_v, Vec3::new(0.0, -10.0, 0.0));
        // assert_approx_eq!(camera.viewport_u, Vec3::new(17.7777777777, 0.0, 0.0));
        assert_approx_eq!(camera.viewport_upper_left, Vec3::new(-8.88888888, 5.0, -5.0));
        assert_approx_eq!(camera.pixel_delta_u, Vec3::new(0.017361111, 0.0, 0.0));
        assert_approx_eq!(camera.pixel_delta_v, Vec3::new(0.0, -0.017361111, 0.0));
        assert_approx_eq!(camera.pixel00_loc, Vec3::new(-8.88020833333, 4.9913194444, -5.0));
    }
}