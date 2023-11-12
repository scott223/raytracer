use crate::linalg::Vec3;
use std::fmt;

#[derive(Debug, Clone, Copy)]
pub struct Ray {
    pub origin: Vec3,
    pub direction: Vec3,
    pub inv_direction: Vec3,

    /// Sign of the X direction. 0 means positive, 1 means negative.
    /// Cached for use in [`AABB`] intersections.
    pub sign_x: usize,

    /// Sign of the Y direction. 0 means positive, 1 means negative.
    /// Cached for use in [`AABB`] intersections.
    pub sign_y: usize,

    /// Sign of the Z direction. 0 means positive, 1 means negative.
    /// Cached for use in [`AABB`] intersections.
    pub sign_z: usize,
}

// Creating a new Ray with a origin and a direction (normalized)
impl Ray {
    pub fn new(origin: Vec3, direction: Vec3) -> Self {
        let dir_norm = direction.normalized();

        let inv_direction = Vec3::new(1.0 / dir_norm.x(), 1.0 / dir_norm.y(), 1.0 / dir_norm.z());

        Ray {
            origin,
            direction: dir_norm,
            inv_direction,
            sign_x: (dir_norm.x() < 0.0) as usize,
            sign_y: (dir_norm.y() < 0.0) as usize,
            sign_z: (dir_norm.z() < 0.0) as usize,
        }
    }

    // finds the location of a point on a ray for a given t
    pub fn at(&self, t: f64) -> Vec3 {
        self.origin + self.direction * t
    }
}

// diplay trait for formatting (e.g. println! or log!)
impl fmt::Display for Ray {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "origin: [{}], direction: [{}]",
            self.origin, self.direction
        )
    }
}

#[cfg(test)]
mod tests {
    use crate::linalg::Vec3;
    use crate::render::ray::Ray;
    use assert_approx_eq::assert_approx_eq;

    #[test_log::test]
    fn test_create_ray() {
        let p = Vec3::new(0.1, 0.2, 0.3);
        let q = Vec3::new(0.2, 0.3, 0.4);

        let r = Ray::new(p, q);

        assert_approx_eq!(r.origin, p);
        assert_approx_eq!(r.direction, q.normalized());
    }

    #[test_log::test]
    fn test_ray_at() {
        let p = Vec3::new(0.0, 0.0, 0.0);
        let q = Vec3::new(0.0, 0.0, 1.0);

        let r = Ray::new(p, q);

        assert_approx_eq!(r.at(5.0), Vec3::new(0.0, 0.0, 5.0));
    }
}
