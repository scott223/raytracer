use crate::vec3::Vec3;

// struct for Probability Density Functions 
#[derive(Debug, Clone, Copy)]
pub struct Pdf {
 //   pub value: Vec3,
}

pub trait Pdf_trait {
    fn generate() -> Vec3;
}

pub struct SpherePDF {
    pub value: f64,
}

impl Pdf_trait for SpherePDF {
    fn generate() -> Vec3 {
        todo!();
    }
}