use serde::{Deserialize, Serialize};
use std::ops::{Add, AddAssign, Div, Mul, MulAssign};

use super::vec3::Vec3;

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct Vec4 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub w: f64,
}

impl Vec4 {
    /// create a new Vec3 with x, y, z, w axis
    pub fn new(x: f64, y: f64, z: f64, w: f64) -> Self {
        Vec4 { x, y, z, w }
    }

    pub fn new_from_vec3(v: Vec3, w: f64) -> Self {
        Vec4 {
            x: v.x(),
            y: v.y(),
            z: v.z(),
            w,
        }
    }

    pub fn to_vec3(self) -> Vec3 {
        Vec3::new(self.x, self.y, self.z)
    }
}

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct Mat4 {
    pub data: [[f64; 4]; 4],
}

impl Mat4 {
    /// Creates a new identity 4x4 matrix
    #[inline(always)]
    #[allow(dead_code)]
    pub fn identity() -> Self {
        Mat4 {
            data: [
                [1., 0., 0., 0.],
                [0., 1., 0., 0.],
                [0., 0., 1., 0.],
                [0., 0., 0., 1.],
            ],
        }
    }

    /// Computes the transpose matrix
    #[inline(always)]
    pub fn transpose(d_x: f64, d_y: f64, d_z: f64) -> Mat4 {
        Mat4 {
            data: [
                [1., 0., 0., d_x],
                [0., 1., 0., d_y],
                [0., 0., 1., d_z],
                [0., 0., 0., 1.],
            ],
        }
    }

    #[inline(always)]
    pub fn scale(s_x: f64, s_y: f64, s_z: f64) -> Mat4 {
        Mat4 {
            data: [
                [s_x, 0., 0., 0.],
                [0., s_y, 0., 0.],
                [0., 0., s_z, 0.],
                [0., 0., 0., 0.],
            ],
        }
    }

    /// Computes the rotation matrix (around the origin)
    /// source: https://www.brainvoyager.com/bv/doc/UsersGuide/CoordsAndTransforms/SpatialTransformationMatrices.html
    #[inline(always)]
    pub fn rotate_x(theta_x: f64) -> Mat4 {
        Mat4 {
            data: [
                [1., 0., 0., 0.],
                [0., theta_x.cos(), theta_x.sin(), 0.],
                [0., -theta_x.sin(), theta_x.cos(), 0.],
                [0., 0., 0., 1.],
            ],
        }
    }

    /// Computes the rotation matrix (around the origin)
    /// source: https://www.brainvoyager.com/bv/doc/UsersGuide/CoordsAndTransforms/SpatialTransformationMatrices.html
    #[inline(always)]
    pub fn rotate_y(theta_y: f64) -> Mat4 {
        Mat4 {
            data: [
                [theta_y.cos(), 0., -theta_y.sin(), 0.],
                [0., 1., 0., 0.],
                [theta_y.sin(), 0., theta_y.cos(), 0.],
                [0., 0., 0., 1.],
            ],
        }
    }

    /// Computes the rotation matrix (around the origin)
    /// source: https://www.brainvoyager.com/bv/doc/UsersGuide/CoordsAndTransforms/SpatialTransformationMatrices.html
    #[inline(always)]
    pub fn rotate_z(theta_z: f64) -> Mat4 {
        Mat4 {
            data: [
                [theta_z.cos(), -theta_z.sin(), 0., 0.],
                [theta_z.sin(), theta_z.cos(), 0., 0.],
                [0., 0., 1., 0.],
                [0., 0., 0., 1.],
            ],
        }
    }

    /// Computes the determinant of a 4x4 matrix
    /// Source: https://docs.piston.rs/graphics/src/vecmath/lib.rs.html
    #[inline(always)]
    #[allow(dead_code)]
    pub fn determinant(&self) -> f64 {
        let mat = self.data;

        mat[0][0] * mat[1][1] * mat[2][2] * mat[3][3]
            + mat[0][0] * mat[1][2] * mat[2][3] * mat[3][1]
            + mat[0][0] * mat[1][3] * mat[2][1] * mat[3][2]
            + mat[0][1] * mat[1][0] * mat[2][3] * mat[3][2]
            + mat[0][1] * mat[1][2] * mat[2][0] * mat[3][3]
            + mat[0][1] * mat[1][3] * mat[2][2] * mat[3][0]
            + mat[0][2] * mat[1][0] * mat[2][1] * mat[3][3]
            + mat[0][2] * mat[1][1] * mat[2][3] * mat[3][0]
            + mat[0][2] * mat[1][3] * mat[2][0] * mat[3][1]
            + mat[0][3] * mat[1][0] * mat[2][2] * mat[3][1]
            + mat[0][3] * mat[1][1] * mat[2][0] * mat[3][2]
            + mat[0][3] * mat[1][2] * mat[2][1] * mat[3][0]
            - mat[0][0] * mat[1][1] * mat[2][3] * mat[3][2]
            - mat[0][0] * mat[1][2] * mat[2][1] * mat[3][3]
            - mat[0][0] * mat[1][3] * mat[2][2] * mat[3][1]
            - mat[0][1] * mat[1][0] * mat[2][2] * mat[3][3]
            - mat[0][1] * mat[1][2] * mat[2][3] * mat[3][0]
            - mat[0][1] * mat[1][3] * mat[2][0] * mat[3][2]
            - mat[0][2] * mat[1][0] * mat[2][3] * mat[3][1]
            - mat[0][2] * mat[1][1] * mat[2][0] * mat[3][3]
            - mat[0][2] * mat[1][3] * mat[2][1] * mat[3][0]
            - mat[0][3] * mat[1][0] * mat[2][1] * mat[3][2]
            - mat[0][3] * mat[1][1] * mat[2][2] * mat[3][0]
            - mat[0][3] * mat[1][2] * mat[2][0] * mat[3][1]
    }

    /// Computes the inverse determinant of a 4x4 matrix
    /// Source: https://docs.piston.rs/graphics/src/vecmath/lib.rs.html
    #[inline(always)]
    #[allow(dead_code)]
    fn inv_determinant(&self) -> f64 {
        1.0 / self.determinant()
    }

    /// Computes the inverse of a 4x4 matrix.
    /// Source: https://docs.piston.rs/graphics/src/vecmath/lib.rs.html
    #[inline(always)]
    #[allow(dead_code)]
    pub fn inverse(&self) -> Mat4 {
        let inv_det: f64 = self.inv_determinant();
        let mat: [[f64; 4]; 4] = self.data;

        Mat4 {
            data: [
                [
                    (mat[1][1] * mat[2][2] * mat[3][3]
                        + mat[1][2] * mat[2][3] * mat[3][1]
                        + mat[1][3] * mat[2][1] * mat[3][2]
                        - mat[1][1] * mat[2][3] * mat[3][2]
                        - mat[1][2] * mat[2][1] * mat[3][3]
                        - mat[1][3] * mat[2][2] * mat[3][1])
                        * inv_det,
                    (mat[0][1] * mat[2][3] * mat[3][2]
                        + mat[0][2] * mat[2][1] * mat[3][3]
                        + mat[0][3] * mat[2][2] * mat[3][1]
                        - mat[0][1] * mat[2][2] * mat[3][3]
                        - mat[0][2] * mat[2][3] * mat[3][1]
                        - mat[0][3] * mat[2][1] * mat[3][2])
                        * inv_det,
                    (mat[0][1] * mat[1][2] * mat[3][3]
                        + mat[0][2] * mat[1][3] * mat[3][1]
                        + mat[0][3] * mat[1][1] * mat[3][2]
                        - mat[0][1] * mat[1][3] * mat[3][2]
                        - mat[0][2] * mat[1][1] * mat[3][3]
                        - mat[0][3] * mat[1][2] * mat[3][1])
                        * inv_det,
                    (mat[0][1] * mat[1][3] * mat[2][2]
                        + mat[0][2] * mat[1][1] * mat[2][3]
                        + mat[0][3] * mat[1][2] * mat[2][1]
                        - mat[0][1] * mat[1][2] * mat[2][3]
                        - mat[0][2] * mat[1][3] * mat[2][1]
                        - mat[0][3] * mat[1][1] * mat[2][2])
                        * inv_det,
                ],
                [
                    (mat[1][0] * mat[2][3] * mat[3][2]
                        + mat[1][2] * mat[2][0] * mat[3][3]
                        + mat[1][3] * mat[2][2] * mat[3][0]
                        - mat[1][0] * mat[2][2] * mat[3][3]
                        - mat[1][2] * mat[2][3] * mat[3][0]
                        - mat[1][3] * mat[2][0] * mat[3][2])
                        * inv_det,
                    (mat[0][0] * mat[2][2] * mat[3][3]
                        + mat[0][2] * mat[2][3] * mat[3][0]
                        + mat[0][3] * mat[2][0] * mat[3][2]
                        - mat[0][0] * mat[2][3] * mat[3][2]
                        - mat[0][2] * mat[2][0] * mat[3][3]
                        - mat[0][3] * mat[2][2] * mat[3][0])
                        * inv_det,
                    (mat[0][0] * mat[1][3] * mat[3][2]
                        + mat[0][2] * mat[1][0] * mat[3][3]
                        + mat[0][3] * mat[1][2] * mat[3][0]
                        - mat[0][0] * mat[1][2] * mat[3][3]
                        - mat[0][2] * mat[1][3] * mat[3][0]
                        - mat[0][3] * mat[1][0] * mat[3][2])
                        * inv_det,
                    (mat[0][0] * mat[1][2] * mat[2][3]
                        + mat[0][2] * mat[1][3] * mat[2][0]
                        + mat[0][3] * mat[1][0] * mat[2][2]
                        - mat[0][0] * mat[1][3] * mat[2][2]
                        - mat[0][2] * mat[1][0] * mat[2][3]
                        - mat[0][3] * mat[1][2] * mat[2][0])
                        * inv_det,
                ],
                [
                    (mat[1][0] * mat[2][1] * mat[3][3]
                        + mat[1][1] * mat[2][3] * mat[3][0]
                        + mat[1][3] * mat[2][0] * mat[3][1]
                        - mat[1][0] * mat[2][3] * mat[3][1]
                        - mat[1][1] * mat[2][0] * mat[3][3]
                        - mat[1][3] * mat[2][1] * mat[3][0])
                        * inv_det,
                    (mat[0][0] * mat[2][3] * mat[3][1]
                        + mat[0][1] * mat[2][0] * mat[3][3]
                        + mat[0][3] * mat[2][1] * mat[3][0]
                        - mat[0][0] * mat[2][1] * mat[3][3]
                        - mat[0][1] * mat[2][3] * mat[3][0]
                        - mat[0][3] * mat[2][0] * mat[3][1])
                        * inv_det,
                    (mat[0][0] * mat[1][1] * mat[3][3]
                        + mat[0][1] * mat[1][3] * mat[3][0]
                        + mat[0][3] * mat[1][0] * mat[3][1]
                        - mat[0][0] * mat[1][3] * mat[3][1]
                        - mat[0][1] * mat[1][0] * mat[3][3]
                        - mat[0][3] * mat[1][1] * mat[3][0])
                        * inv_det,
                    (mat[0][0] * mat[1][3] * mat[2][1]
                        + mat[0][1] * mat[1][0] * mat[2][3]
                        + mat[0][3] * mat[1][1] * mat[2][0]
                        - mat[0][0] * mat[1][1] * mat[2][3]
                        - mat[0][1] * mat[1][3] * mat[2][0]
                        - mat[0][3] * mat[1][0] * mat[2][1])
                        * inv_det,
                ],
                [
                    (mat[1][0] * mat[2][2] * mat[3][1]
                        + mat[1][1] * mat[2][0] * mat[3][2]
                        + mat[1][2] * mat[2][1] * mat[3][0]
                        - mat[1][0] * mat[2][1] * mat[3][2]
                        - mat[1][1] * mat[2][2] * mat[3][0]
                        - mat[1][2] * mat[2][0] * mat[3][1])
                        * inv_det,
                    (mat[0][0] * mat[2][1] * mat[3][2]
                        + mat[0][1] * mat[2][2] * mat[3][0]
                        + mat[0][2] * mat[2][0] * mat[3][1]
                        - mat[0][0] * mat[2][2] * mat[3][1]
                        - mat[0][1] * mat[2][0] * mat[3][2]
                        - mat[0][2] * mat[2][1] * mat[3][0])
                        * inv_det,
                    (mat[0][0] * mat[1][2] * mat[3][1]
                        + mat[0][1] * mat[1][0] * mat[3][2]
                        + mat[0][2] * mat[1][1] * mat[3][0]
                        - mat[0][0] * mat[1][1] * mat[3][2]
                        - mat[0][1] * mat[1][2] * mat[3][0]
                        - mat[0][2] * mat[1][0] * mat[3][1])
                        * inv_det,
                    (mat[0][0] * mat[1][1] * mat[2][2]
                        + mat[0][1] * mat[1][2] * mat[2][0]
                        + mat[0][2] * mat[1][0] * mat[2][1]
                        - mat[0][0] * mat[1][2] * mat[2][1]
                        - mat[0][1] * mat[1][0] * mat[2][2]
                        - mat[0][2] * mat[1][1] * mat[2][0])
                        * inv_det,
                ],
            ],
        }
    }
}

/// Multiply Mat4 with a Vec3
impl Mul<Vec4> for Mat4 {
    type Output = Vec4;

    fn mul(self, v: Vec4) -> Vec4 {
        Vec4::new(
            self.data[0][0] * v.x
                + self.data[0][1] * v.y
                + self.data[0][2] * v.z
                + self.data[0][3] * v.w,
            self.data[1][0] * v.x
                + self.data[1][1] * v.y
                + self.data[1][2] * v.z
                + self.data[1][3] * v.w,
            self.data[2][0] * v.x
                + self.data[2][1] * v.y
                + self.data[2][2] * v.z
                + self.data[2][3] * v.w,
            self.data[3][0] * v.x
                + self.data[3][1] * v.y
                + self.data[3][2] * v.z
                + self.data[3][3] * v.w,
        )
    }
}

/// Add Mat4 to a Mat4
impl Add<Mat4> for Mat4 {
    type Output = Mat4;

    fn add(self, m: Mat4) -> Mat4 {
        Mat4 {
            data: [
                [
                    self.data[0][0] + m.data[0][0],
                    self.data[0][1] + m.data[0][1],
                    self.data[0][2] + m.data[0][2],
                    self.data[0][3] + m.data[0][3],
                ],
                [
                    self.data[1][0] + m.data[1][0],
                    self.data[1][1] + m.data[1][1],
                    self.data[1][2] + m.data[1][2],
                    self.data[1][3] + m.data[1][3],
                ],
                [
                    self.data[2][0] + m.data[2][0],
                    self.data[2][1] + m.data[2][1],
                    self.data[2][2] + m.data[2][2],
                    self.data[2][3] + m.data[2][3],
                ],
                [
                    self.data[3][0] + m.data[3][0],
                    self.data[3][1] + m.data[3][1],
                    self.data[3][2] + m.data[3][2],
                    self.data[3][3] + m.data[3][3],
                ],
            ],
        }
    }
}

/// Add assign Mat4 to self
impl AddAssign<Mat4> for Mat4 {
    fn add_assign(&mut self, m: Mat4) {
        self.data = [
            [
                self.data[0][0] + m.data[0][0],
                self.data[0][1] + m.data[0][1],
                self.data[0][2] + m.data[0][2],
                self.data[0][3] + m.data[0][3],
            ],
            [
                self.data[1][0] + m.data[1][0],
                self.data[1][1] + m.data[1][1],
                self.data[1][2] + m.data[1][2],
                self.data[1][3] + m.data[1][3],
            ],
            [
                self.data[2][0] + m.data[2][0],
                self.data[2][1] + m.data[2][1],
                self.data[2][2] + m.data[2][2],
                self.data[2][3] + m.data[2][3],
            ],
            [
                self.data[3][0] + m.data[3][0],
                self.data[3][1] + m.data[3][1],
                self.data[3][2] + m.data[3][2],
                self.data[3][3] + m.data[3][3],
            ],
        ];
    }
}

/// mul assign Mat4 to self
impl MulAssign<Mat4> for Mat4 {
    fn mul_assign(&mut self, m: Mat4) {
        self.data = [
            [
                self.data[0][0] * m.data[0][0],
                self.data[0][1] * m.data[0][1],
                self.data[0][2] * m.data[0][2],
                self.data[0][3] * m.data[0][3],
            ],
            [
                self.data[1][0] * m.data[1][0],
                self.data[1][1] * m.data[1][1],
                self.data[1][2] * m.data[1][2],
                self.data[1][3] * m.data[1][3],
            ],
            [
                self.data[2][0] * m.data[2][0],
                self.data[2][1] * m.data[2][1],
                self.data[2][2] * m.data[2][2],
                self.data[2][3] * m.data[2][3],
            ],
            [
                self.data[3][0] * m.data[3][0],
                self.data[3][1] * m.data[3][1],
                self.data[3][2] * m.data[3][2],
                self.data[3][3] * m.data[3][3],
            ],
        ];
    }
}

// Divide by a f64
impl Div<f64> for Mat4 {
    type Output = Mat4;

    fn div(self, q: f64) -> Mat4 {
        Mat4 {
            data: [
                [
                    self.data[0][0] / q,
                    self.data[0][1] / q,
                    self.data[0][2] / q,
                    self.data[0][3] / q,
                ],
                [
                    self.data[1][0] / q,
                    self.data[1][1] / q,
                    self.data[1][2] / q,
                    self.data[1][3] / q,
                ],
                [
                    self.data[2][0] / q,
                    self.data[2][1] / q,
                    self.data[2][2] / q,
                    self.data[2][3] / q,
                ],
                [
                    self.data[3][0] / q,
                    self.data[3][1] / q,
                    self.data[3][2] / q,
                    self.data[3][3] / q,
                ],
            ],
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::linalg::{
        mat4::{Mat4, Vec4},
        vec3::Vec3,
    };
    use assert_approx_eq::assert_approx_eq;

    // Test creating a new vector with three rows, taking floats as an argument
    #[test_log::test]
    fn test_new_identity() {
        let m = Mat4::identity();

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
        let m: Mat4 = Mat4::identity();
        let v: Vec3 = Vec3::new(1., 2., 3.);
        let vv: Vec4 = Vec4::new_from_vec3(v, 1.);

        let p: Vec4 = m * vv;

        assert_approx_eq!(p.x, 1.0);
        assert_approx_eq!(p.y, 2.0);
        assert_approx_eq!(p.z, 3.0);
        assert_approx_eq!(p.w, 1.0);

        let mm: Mat4 = Mat4 {
            data: [
                [1., 2., 3., 4.],
                [2., 3., 4., 5.],
                [3., 4., 5., 6.],
                [0., 0., 0., 1.],
            ],
        };

        let pp = mm * vv;

        assert_approx_eq!(pp.x, 1.0 * 1.0 + 2. * 2.0 + 3. * 3. + 4. * 1.);
        assert_approx_eq!(pp.y, 2.0 * 1.0 + 3. * 2.0 + 4. * 3. + 5. * 1.);
        assert_approx_eq!(pp.z, 3.0 * 1.0 + 4. * 2.0 + 5. * 3. + 6. * 1.);
        assert_approx_eq!(pp.w, 0.0 * 1.0 + 0. * 2.0 + 0. * 3. + 1. * 1.);
    }

    #[test_log::test]
    fn test_divide_f64() {
        let m: Mat4 = Mat4::identity();
        let f: f64 = 2.0;
        let p: Mat4 = m / f;

        assert_approx_eq!(p.data[0][0], 0.5);
        assert_approx_eq!(p.data[0][1], 0.);
        assert_approx_eq!(p.data[0][2], 0.);
        assert_approx_eq!(p.data[0][3], 0.);

        assert_approx_eq!(p.data[1][0], 0.);
        assert_approx_eq!(p.data[1][1], 0.5);
        assert_approx_eq!(p.data[1][2], 0.);
        assert_approx_eq!(p.data[1][3], 0.);

        assert_approx_eq!(p.data[2][0], 0.);
        assert_approx_eq!(p.data[2][1], 0.);
        assert_approx_eq!(p.data[2][2], 0.5);
        assert_approx_eq!(p.data[2][3], 0.);

        assert_approx_eq!(p.data[3][0], 0.);
        assert_approx_eq!(p.data[3][1], 0.);
        assert_approx_eq!(p.data[3][2], 0.);
        assert_approx_eq!(p.data[3][3], 0.5);
    }

    #[test_log::test]
    fn test_add() {
        let m: Mat4 = Mat4::identity();
        let v: Mat4 = Mat4::identity();

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
            data: [
                [1., 2., 3., 4.],
                [2., 3., 4., 5.],
                [3., 4., 5., 6.],
                [0., 0., 0., 1.],
            ],
        };

        let vv: Mat4 = Mat4 {
            data: [
                [1., 2., 3., 4.],
                [2., 3., 4., 5.],
                [3., 4., 5., 6.],
                [0., 0., 0., 1.],
            ],
        };

        let pp: Mat4 = mm + vv;

        assert_approx_eq!(pp.data[0][0], 2.0);
        assert_approx_eq!(pp.data[1][1], 6.0);
        assert_approx_eq!(pp.data[2][2], 10.0);
        assert_approx_eq!(pp.data[3][3], 2.0);
    }

    #[test_log::test]
    fn test_determinant() {
        let a: Mat4 = Mat4 {
            data: [
                [1., 1., 1., -1.],
                [1., 1., -1., 1.],
                [1., -1., 1., 1.],
                [-1., 1., 1., 1.],
            ],
        };

        let determinant = a.determinant();
        assert_approx_eq!(determinant, -16.0);
    }

    #[test_log::test]
    fn test_inverse() {
        let a: Mat4 = Mat4 {
            data: [
                [1., 1., 1., -1.],
                [1., 1., -1., 1.],
                [1., -1., 1., 1.],
                [-1., 1., 1., 1.],
            ],
        };

        let inv_a = a.inverse();

        assert_approx_eq!(inv_a.data[0][0], 0.25);
        assert_approx_eq!(inv_a.data[0][1], 0.25);
        assert_approx_eq!(inv_a.data[0][2], 0.25);
        assert_approx_eq!(inv_a.data[0][3], -0.25);

        assert_approx_eq!(inv_a.data[1][0], 0.25);
        assert_approx_eq!(inv_a.data[1][1], 0.25);
        assert_approx_eq!(inv_a.data[1][2], -0.25);
        assert_approx_eq!(inv_a.data[1][3], 0.25);

        assert_approx_eq!(inv_a.data[2][0], 0.25);
        assert_approx_eq!(inv_a.data[2][1], -0.25);
        assert_approx_eq!(inv_a.data[2][2], 0.25);
        assert_approx_eq!(inv_a.data[2][3], 0.25);

        assert_approx_eq!(inv_a.data[3][0], -0.25);
        assert_approx_eq!(inv_a.data[3][1], 0.25);
        assert_approx_eq!(inv_a.data[3][2], 0.25);
        assert_approx_eq!(inv_a.data[3][3], 0.25);
    }
}
