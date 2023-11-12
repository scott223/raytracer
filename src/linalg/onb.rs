use super::vec3::Vec3;

// struct for Orthonormal Bases
#[derive(Debug, Clone, Copy)]
pub struct Onb {
    pub axis: [Vec3; 3],
}

impl Onb {
    pub fn u(self) -> Vec3 {
        self.axis[0]
    }

    pub fn v(self) -> Vec3 {
        self.axis[1]
    }

    pub fn w(self) -> Vec3 {
        self.axis[2]
    }

    #[allow(dead_code)]
    pub fn local_abc(self, a: f64, b: f64, c: f64) -> Vec3 {
        self.u() * a + self.v() * b + self.w() * c
    }

    pub fn local_vec(self, vec: Vec3) -> Vec3 {
        self.u() * vec.x() + self.v() * vec.y() + self.w() * vec.z()
    }

    pub fn build_from_w(w: Vec3) -> Self {
        let unit_w = w.normalized();
        let a = if unit_w.x().abs() > 0.9 {
            Vec3::new(0., 1., 0.)
        } else {
            Vec3::new(1.0, 0., 0.)
        };

        let v = unit_w.cross(&a).normalized();
        let u = unit_w.cross(&v);

        Self {
            axis: [u, v, unit_w],
        }
    }
}
