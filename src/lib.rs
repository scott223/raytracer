// Exteral imports
use std::error::Error;

// Specfic imports
mod vec3;
mod ray;
use crate::vec3::Vec3;
use crate::ray::Ray;

pub fn render() -> Result<(), Box<dyn Error>> {
    let o: Vec3 = Vec3::new(0 as f64,0 as f64,0 as f64);
    let pixel_00: Vec3 = Vec3::new(1.1, 1.1, -2.0);
    let ray_00: Ray = Ray::new(o, pixel_00);

    Ok(())
}

#[cfg(test)]
mod tests {
    use std::error::Error;
    use crate::render;
    
    #[test_log::test]
    fn test_render_full_scene() -> Result <(), Box<dyn Error>> {
        render()?;
        Ok(())
    }
}