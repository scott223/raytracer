
use std::ops::{Add, Neg, Sub, Mul, Div};
use serde::{Serialize, Deserialize};

use crate::vec3::Vec3;

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct Vec4 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub w: f64,
}

impl Vec4 {
    pub fn new_from_vec3(v: Vec3, w: f64) -> Self {
        Vec4 {
            x: v.x(),
            y: v.y(),
            z: v.z(),
            w: w,
        }
    }

    pub fn to_vec3(&self) -> Vec3 {
        Vec3::new(
            self.x,
            self.y, 
            self.z,
        )
    }

    pub fn new(x: f64, y: f64, z: f64, w: f64) -> Self {
        Vec4 {
            x,
            y,
            z,
            w,
        }
    }
}

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct Mat4 {
    pub data: [[f64; 4]; 4],
}

impl Mat4 {
    pub fn new_identity() -> Self {
        Mat4 {
            data: [[1., 0., 0., 0.],
                  [0., 1., 0., 0.],
                  [0., 0., 1., 0.],
                  [0., 0., 0., 1.]]
        }
    }
}

// multiply with a Vec3 
impl Mul<Vec4> for Mat4 {
    type Output = Vec4;

    fn mul(self, v: Vec4) -> Vec4 {
        Vec4::new(
            self.data[0][0] * v.x + self.data[0][1] * v.y + self.data[0][2] * v.z + self.data[0][3] * v.w,
            self.data[1][0] * v.x + self.data[1][1] * v.y + self.data[1][2] * v.z + self.data[1][3] * v.w,
            self.data[2][0] * v.x + self.data[2][1] * v.y + self.data[2][2] * v.z + self.data[2][3] * v.w,
            self.data[3][0] * v.x + self.data[3][1] * v.y + self.data[3][2] * v.z + self.data[3][3] * v.w,
        )
    }
}

// multiply with a Vec3 
impl Add<Mat4> for Mat4 {
    type Output = Mat4;

    fn add(self, m: Mat4) -> Mat4 {
        Mat4 {
            data: [[self.data[0][0] + m.data[0][0], self.data[0][1] + m.data[0][1], self.data[0][2] + m.data[0][2], self.data[0][3] + m.data[0][3]],
                  [self.data[1][0] + m.data[1][0], self.data[1][1] + m.data[1][1], self.data[1][2] + m.data[1][2], self.data[1][3] + m.data[1][3]],
                  [self.data[2][0] + m.data[2][0], self.data[2][1] + m.data[2][1], self.data[2][2] + m.data[2][2], self.data[2][3] + m.data[2][3]],
                  [self.data[3][0] + m.data[3][0], self.data[3][1] + m.data[3][1], self.data[3][2] + m.data[3][2], self.data[3][3] + m.data[3][3]]]
        }
    }
}

#[cfg(test)]
mod tests {
    use assert_approx_eq::assert_approx_eq;
    use crate::vec3::Vec3;
    use crate::mat4::Vec4;
    use crate::mat4::Mat4;

    // Test creating a new vector with three rows, taking floats as an argument
    #[test_log::test]
    fn test_new_identity() {
        let m = Mat4::new_identity();

        assert_approx_eq!(m.data[0][0], 1.0);
        assert_approx_eq!(m.data[1][1], 1.0);
        assert_approx_eq!(m.data[2][2], 1.0);
        assert_approx_eq!(m.data[3][3], 1.0);

        assert_approx_eq!(m.data[1][0], 0.);
        assert_approx_eq!(m.data[2][3], 0.);
        assert_approx_eq!(m.data[1][3], 0.);
        assert_approx_eq!(m.data[3][1], 0.);
    }

    #[test_log::test]
    fn test_vec4() {
        let v = Vec3::new(1., 2., 3.);
        let vv = Vec4::new_from_vec3(v, 1.);

        assert_approx_eq!(vv.x, 1.0);
        assert_approx_eq!(vv.y, 2.0);
        assert_approx_eq!(vv.z, 3.0);
        assert_approx_eq!(vv.w, 1.0);

        let p = vv.to_vec3();

        assert_approx_eq!(p.x(), 1.0);
        assert_approx_eq!(p.y(), 2.0);
        assert_approx_eq!(p.z(), 3.0);  
    }

    #[test_log::test]
    fn test_multiply_vec3() {
        let m: Mat4 = Mat4::new_identity();
        let v: Vec3 = Vec3::new(1., 2., 3.);
        let vv: Vec4 = Vec4::new_from_vec3(v, 1.);

        let p: Vec4 = m * vv;

        assert_approx_eq!(p.x, 1.0);
        assert_approx_eq!(p.y, 2.0);
        assert_approx_eq!(p.z, 3.0);
        assert_approx_eq!(p.w, 1.0);

       let mm: Mat4 = Mat4 {
            data: [[1., 2., 3., 4.],
                  [2., 3., 4., 5.],
                  [3., 4., 5., 6.],
                  [0., 0., 0., 1.]]
        };

        let pp = mm * vv;

        assert_approx_eq!(pp.x, 1.0 * 1.0 + 2. * 2.0 + 3. * 3. + 4. * 1.);
        assert_approx_eq!(pp.y, 2.0 * 1.0 + 3. * 2.0 + 4. * 3. + 5. * 1.);
        assert_approx_eq!(pp.z, 3.0 * 1.0 + 4. * 2.0 + 5. * 3. + 6. * 1.);
        assert_approx_eq!(pp.w, 0.0 * 1.0 + 0. * 2.0 + 0. * 3. + 1. * 1.);
    }

    #[test_log::test]
    fn test_add() {

        let m: Mat4 = Mat4::new_identity();
        let v: Mat4 = Mat4::new_identity();

        let p = m + v;

        assert_approx_eq!(p.data[0][0], 2.0);
        assert_approx_eq!(p.data[1][1], 2.0);
        assert_approx_eq!(p.data[2][2], 2.0);
        assert_approx_eq!(p.data[3][3], 2.0);

        assert_approx_eq!(p.data[1][0], 0.);
        assert_approx_eq!(p.data[2][3], 0.);
        assert_approx_eq!(p.data[1][3], 0.);
        assert_approx_eq!(p.data[3][1], 0.);

        let mm: Mat4 = Mat4 {
            data: [[1., 2., 3., 4.],
                  [2., 3., 4., 5.],
                  [3., 4., 5., 6.],
                  [0., 0., 0., 1.]]
        };

        let vv: Mat4 = Mat4 {
            data: [[1., 2., 3., 4.],
                  [2., 3., 4., 5.],
                  [3., 4., 5., 6.],
                  [0., 0., 0., 1.]]
        };

        let pp: Mat4 = mm + vv;

        assert_approx_eq!(pp.data[0][0], 2.0);
        assert_approx_eq!(pp.data[1][1], 6.0);
        assert_approx_eq!(pp.data[2][2], 10.0);
        assert_approx_eq!(pp.data[3][3], 2.0);

    }
}