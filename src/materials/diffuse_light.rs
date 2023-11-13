use serde::{Deserialize, Serialize};

use crate::{elements::HitRecord, render::Color, render::Ray};

use super::Emmits;

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct DiffuseLight {
    // albedo is here the amount of light emitted (can be more than >1.)
    pub albedo: Color,
}

// create a new Diffuse light material
impl DiffuseLight {
    pub fn new(albedo: Color) -> DiffuseLight {
        DiffuseLight { albedo }
    }
}

impl Emmits for DiffuseLight {
    fn emitted(&self, _ray: &Ray, hit_record: &HitRecord) -> Option<Color> {
        //only return the color if we hit the front side of the emmitting element
        if hit_record.front_face {
            Some(self.albedo)
        } else {
            None
        }
    }
}
