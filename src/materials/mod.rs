mod diffuse_light;
pub use diffuse_light::DiffuseLight;

mod lambertian;
pub use lambertian::Lambertian;

mod metal;
pub use metal::Metal;

mod dielectric;
pub use dielectric::Dielectric;

use rand::Rng;
use serde::{Deserialize, Serialize};

use crate::{
    elements::HitRecord,
    render::Color,
    render::{Pdf, Ray},
};

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub enum Material {
    Lambertian(Lambertian),
    Metal(Metal),
    Dielectric(Dielectric),
    DiffuseLight(DiffuseLight),
}

//need to specify lifetime for the pdf, as there might be a mixed pdf in there with a reference to other pdfs
pub struct ScatterRecord<'a> {
    pub attenuation: Color,
    pub pdf: Pdf<'a>,
}

// trait for a material that scatters
pub trait Scatterable {
    fn scatter(
        &self,
        ray: &Ray,
        hit_record: &HitRecord,
        rng: &mut impl Rng,
    ) -> Option<ScatterRecord>;
}

// link the trait implementation to the materials
impl Scatterable for Material {
    fn scatter(
        &self,
        ray: &Ray,
        hit_record: &HitRecord,
        rng: &mut impl Rng,
    ) -> Option<ScatterRecord> {
        match self {
            Material::Lambertian(l) => l.scatter(ray, hit_record, rng),
            _ => None,
        }
    }
}

// reflectance of a material
pub struct ReflectRecord {
    pub attenuation: Color,
    // note that Ray is an option here, if the material absorbs the ray due to the fuzz factor
    pub ray: Option<Ray>,
}

// trait for a material that reflects
pub trait Reflects {
    fn reflect(
        &self,
        ray: &Ray,
        hit_record: &HitRecord,
        rng: &mut impl Rng,
    ) -> Option<ReflectRecord>;
}

// link the trait implementation to the materials
impl Reflects for Material {
    fn reflect(
        &self,
        ray: &Ray,
        hit_record: &HitRecord,
        rng: &mut impl Rng,
    ) -> Option<ReflectRecord> {
        match self {
            Material::Metal(m) => m.reflect(ray, hit_record, rng),
            _ => None,
        }
    }
}

pub struct RefractRecord {
    pub attenuation: Color,
    pub ray: Ray,
}

// trait for a material that refracts
pub trait Refracts {
    fn refract(
        &self,
        ray: &Ray,
        hit_record: &HitRecord,
        rng: &mut impl Rng,
    ) -> Option<RefractRecord>;
}

// link the trait implementation to the materials
impl Refracts for Material {
    fn refract(
        &self,
        ray: &Ray,
        hit_record: &HitRecord,
        rng: &mut impl Rng,
    ) -> Option<RefractRecord> {
        match self {
            Material::Dielectric(d) => d.refract(ray, hit_record, rng),
            _ => None,
        }
    }
}

//trait for a material that emmits, will just return a color and nothing fancy
pub trait Emmits {
    fn emitted(&self, ray: &Ray, hit_record: &HitRecord) -> Option<Color>;
}

impl Emmits for Material {
    fn emitted(&self, ray: &Ray, hit_record: &HitRecord) -> Option<Color> {
        match self {
            Material::DiffuseLight(d) => d.emitted(ray, hit_record),
            _ => None,
        }
    }
}
