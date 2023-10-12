use std::ops::{Add, AddAssign, Div, Mul, Sub};

use serde::{Deserialize, Serialize};

pub struct Rgb {
    pub r: u8,
    pub g: u8,
    pub b: u8,
}

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct Color {
    pub r: f64,
    pub g: f64,
    pub b: f64,
}

impl Color {
    // creates a new color
    pub fn new(r: f64, g: f64, b: f64) -> Self {
        Color { 
            r,
            g,
            b,
        }
    }

    // divides the color by the number of samples
    pub fn divide_by_samples(&self, samples: usize) -> Self {
        Color::new(
            self.r / samples as f64,
            self.g / samples as f64,
            self.b / samples as f64,
        )
    }

    // returns its own value with max 1.0 and min 0.0
    pub fn clamp(&self) -> Self {
        Color::new(
            self.r.min(1.0).max(0.0),
            self.g.min(1.0).max(0.0),
            self.b.min(1.0).max(0.0),
        )
    }

    //apply gamma correction
    pub fn linear_to_gamma(&self, correction: f64) -> Self {
        Color::new(self.r.powf(1.0/correction), self.g.powf(1.0/correction), self.b.powf(1.0/correction))
    }

    //to Rgb, first clamp to max 1.0 and min 0.0. so we dont overflow
    pub fn to_rgb(&self) -> Rgb {
        let c = self.clamp();
        Rgb {
            r: (c.r * 255.0) as u8,
            g: (c.g * 255.0) as u8,
            b: (c.b * 255.0) as u8,
        }
    }
}

// + operator
impl Add for Color {
    type Output = Color;

    fn add(self, w: Color) -> Color {
        Color {
            r: self.r + w.r,
            g: self.g + w.g,
            b: self.b + w.b,
        }
    }
}

// + operator
impl AddAssign for Color {
    fn add_assign(&mut self, w: Color) {
        *self = Self {
            r: self.r + w.r,
            g: self.g + w.g,
            b: self.b + w.b,
        }
    }
}

// - operator
impl Sub for Color {
    type Output = Color;

    fn sub(self, w: Color) -> Color {
        Color {
            r: self.r - w.r,
            g: self.g - w.g,
            b: self.b - w.b,
        }
    }
}

// multiply with a color
impl Mul for Color {
    type Output = Color;

    fn mul(self, q: Color) -> Color {
        Color {
            r: self.r * q.r,
            g: self.g * q.g,
            b: self.b * q.b,
        }
    }
}

// multiply with a f64
impl Mul<f64> for Color {
    type Output = Color;

    fn mul(self, q: f64) -> Color {
        Color {
            r: self.r * q,
            g: self.g * q,
            b: self.b * q,
        }
    }
}

// divide by a u8
impl Div<u8> for Color {
    type Output = Color;

    fn div(self, q: u8) -> Color {
        Color {
            r: self.r / q as f64,
            g: self.g / q as f64,
            b: self.b / q as f64,
        }
    }
}

// divide by a u8
impl Div<f64> for Color {
    type Output = Color;

    fn div(self, q: f64) -> Color {
        Color {
            r: self.r / q as f64,
            g: self.g / q as f64,
            b: self.b / q as f64,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::color::Color;
    use assert_approx_eq::assert_approx_eq;

    #[test_log::test]
    fn test_create_color() {
        let c = Color::new(0.5, 0.4, 0.3);
        assert_approx_eq!(c.r, 0.5);
        assert_approx_eq!(c.g, 0.4);
        assert_approx_eq!(c.b, 0.3);
    }

    #[test_log::test]
    fn test_add() {
        let p = Color::new(0.1, 0.2, 0.3);
        let q = Color::new(0.2, 0.3, 0.4);
        let r = p + q;
        assert_approx_eq!(r.r, 0.3);
        assert_approx_eq!(r.g, 0.5);
        assert_approx_eq!(r.b, 0.7);
    }

    #[test_log::test]
    fn test_sub() {
        let p = Color::new(0.1, 0.2, 0.3);
        let q = Color::new(0.2, 0.3, 0.4);
        let r = p - q;
        assert_approx_eq!(r.r, -0.1);
        assert_approx_eq!(r.g, -0.1);
        assert_approx_eq!(r.b, -0.1);
    }

    #[test_log::test]
    fn test_div() {
        let p = Color::new(0.4, 0.2, 0.8);
        let q: u8 = 2;
        let r = p / q;
        assert_approx_eq!(r.r, 0.2);
        assert_approx_eq!(r.g, 0.1);
        assert_approx_eq!(r.b, 0.4);
    }

    #[test_log::test]
    fn test_clamp_color() {
        let c = Color::new(0.5, 0.4, -2.3);
        let d = Color::new(0.1, 1.0, 0.1);
        let e = c + d;
        println!("e: {:?}",e);
        let f = e.clamp();

        assert_approx_eq!(f.r, 0.6);
        assert_approx_eq!(f.g, 1.0);
        assert_approx_eq!(f.b, 0.0);
    }

    #[test_log::test]
    fn test_gamma_color() {
        let correction: f64 = 2.2;
        let c = Color::new(0.5, 0.4, 0.3).linear_to_gamma(correction);

        assert_approx_eq!(c.r, 0.5_f64.powf(1.0/correction));
        assert_approx_eq!(c.g, 0.4_f64.powf(1.0/correction));
        assert_approx_eq!(c.b, 0.3_f64.powf(1.0/correction));
    }

    #[test_log::test]
    fn test_to_rgb() {
        let c = Color::new(1.2, 0.0, 0.5).to_rgb();

        assert_eq!(c.r, 255);
        assert_eq!(c.g, 0);
        assert_eq!(c.b, 127);
    }
}