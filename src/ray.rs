use crate::vec3::Vec3;
use std::fmt;

#[derive(Debug, Clone, Copy)]
pub struct Ray {
    pub origin: Vec3,
    pub direction: Vec3,
}

// Creating a new Ray with a origin and a direction (normalized)
impl Ray {
    pub fn new(origin: Vec3, direction: Vec3) -> Self {
        Ray {
            origin,
            direction: direction.normalized(),
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
    use crate::ray::Ray;
    use crate::vec3::Vec3;
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
