use crate::vec3::Vec3;
use std::fmt;

#[cfg(test)]
use assert_approx_eq::assert_approx_eq;

#[derive(Debug, Clone, Copy)]
pub struct Ray {
    origin: Vec3,
    direction: Vec3,
}

impl Ray {
    // Creating a new Ray with a origin and a direction (normalized)
    pub fn new(o_val: Vec3, d_val: Vec3) -> Self {
        Ray{
            origin: o_val,
            direction: d_val.normalized(),
        }
    }

    pub fn at(&self, t: f64) -> Vec3 {
        self.origin + self.direction * t
    }
}

impl fmt::Display for Ray {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
      write!(f, "origin: [{}], direction: [{}]", self.origin, self.direction)
    }
}

#[cfg(test)]
mod tests {
    use assert_approx_eq::assert_approx_eq;
    use crate::vec3::Vec3;
    use crate::ray::Ray;
    
    #[test_log::test]
    fn test_create_ray() {
        let p = Vec3::new(0.1, 0.2, 0.3);
        let q = Vec3::new(0.2, 0.3, 0.4);

        let r = Ray::new(p, q);

        assert_approx_eq!(r.origin, p);
        assert_approx_eq!(r.direction, q.normalized());
    }
}