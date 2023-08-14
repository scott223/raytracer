#[derive(Debug, Clone, Copy)]
pub struct Vec3 {
    x: f64,
    y: f64,
    z: f64,
}

impl Vec3 {
    // Creating a new vector
    pub fn new(x_val: f64, y_val: f64, z_val: f64) -> Self {
        Vec3 {
            x: x_val,
            y: y_val,
            z: z_val,
        }
    }

    // Exposing the x, y, z coordinates through functions
    pub fn x(&self) -> f64 {
        self.x
    }

    pub fn y(&self) -> f64 {
        self.y
    }

    pub fn z(&self) -> f64 {
        self.z
    }

    pub fn distance(&self, w: &Vec3) -> f64 {
        let dx = f64::powf(self.x() - w.x(), 2.0);
        let dy = f64::powf(self.y() - w.y(), 2.0);
        let dz = f64::powf(self.z() - w.z(), 2.0);
        (dx + dy + dz).sqrt()
    }

    //pub fn normalize(&self) -> Self {

    //}
}

#[cfg(test)]
mod tests {
    use assert_approx_eq::assert_approx_eq;
    use crate::vec3::Vec3;

    // Check creating a new vector with three rows, taking floats as an argument
    #[test_log::test]
    fn test_vec3_create() {
        let p: Vec3 = Vec3::new(0.1, 0.2, 0.3);
        assert_approx_eq!(p.x(), 0.1);
        assert_approx_eq!(p.y(), 0.2);
        assert_approx_eq!(p.z(), 0.3);

        let q = Vec3::new(0.2, 0.3, 0.4);
        assert_approx_eq!(q.x(), 0.2);
        assert_approx_eq!(q.y(), 0.3);
        assert_approx_eq!(q.z(), 0.4);
    }
    #[test_log::test]
    fn test_vec3_distance() {
        let p: Vec3 = Vec3::new(1.0, 0.0, 5.0);
        let q: Vec3 = Vec3::new(0.0, 2.0, 4.0);
        assert_approx_eq!(p.distance(&q), 2.44948974278);
    }
}