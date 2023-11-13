use serde::{Deserialize, Serialize};

use super::Vec3;

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
