use std::ops::{Neg, Add, Sub, Mul, Div};
use std::fmt;
use rand::Rng;

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

    pub fn new_random(val_min: f64, val_max: f64) -> Self {
        Vec3 {
            x: rand::thread_rng().gen_range(val_min..val_max),
            y: rand::thread_rng().gen_range(val_min..val_max),
            z: rand::thread_rng().gen_range(val_min..val_max),           
        }
    }

    pub fn new_random_unit_sphere() -> Self {
        loop {
            let rv = Vec3::new_random(-1.0,1.0);
            if rv.length_squared() < 1.0 {
                return rv
            }
        }
    }

    pub fn new_random_unit_vector() -> Self {
        Vec3::new_random_unit_sphere().normalized()
    }

    pub fn new_random_on_hemisphere(normal: &Vec3) -> Self {
        let rv = Vec3::new_random_unit_vector();
        if rv.dot(&normal) > 0.0 { //in the same hemisphere
            return rv 
        } else {
            return -rv;
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

    // fn distance
    // calculate the Euclidean distance between two points
    // d(p, w) = sqrt((px - wx)^2 + (py - wy)^2 + (pz - wz)^2)
    pub fn distance(&self, w: &Vec3) -> f64 {
        let dx = f64::powf(self.x - w.x(), 2.0);
        let dy = f64::powf(self.y - w.y(), 2.0);
        let dz = f64::powf(self.z - w.z(), 2.0);
        (dx + dy + dz).sqrt()
    }

    // fn length
    // calculate the vector (length) between point p and the origin (0,0,0)
    pub fn length(&self) -> f64 {
        self.distance(&Vec3::new(0.0, 0.0, 0.0))
    }

    pub fn length_squared(&self) -> f64 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    // fn abs
    // returns absolute value of vector, which is the same as its length (distance from 0.0)
    pub fn abs(&self) -> f64 {
        self.length()
    }

    // vector dot multiplication
    pub fn dot(&self, w: &Vec3) -> f64 {
        self.x * w.x + self.y * w.y + self.z * w.z
    }

    // fn normalized
    // return the normalized vector (= unit vector)
    pub fn normalized(&self) -> Self {
        let l = self.length();
        Vec3::new(self.x()/l, self.y()/l, self.z()/l)
    }

}

// display trait
impl fmt::Display for Vec3 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
      write!(f, "[{:.3}, {:.3}, {:.3}]", self.x, self.y, self.z)
    }
}

// negative value
impl Neg for Vec3 {
    type Output = Vec3;

    fn neg(self) -> Vec3 {
        Vec3 {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}

// + operator
impl Add for Vec3 {
    type Output = Vec3;

    fn add(self, w: Vec3) -> Vec3 {
        Vec3 {
            x: self.x + w.x(),
            y: self.y + w.y(),
            z: self.z + w.z(),
        }
    }
}

// - operator
impl Sub for Vec3 {
    type Output = Vec3;

    fn sub(self, w: Vec3) -> Vec3 {
        Vec3 {
            x: self.x - w.x(),
            y: self.y - w.y(),
            z: self.z - w.z(),
        }
    }
}

// multiply with a Vec3 
impl Mul<Vec3> for Vec3 {
    type Output = Vec3;

    fn mul(self, v: Vec3) -> Vec3 {
        Vec3{
            x: self.x * v.x(),
            y: self.y * v.y(),
            z: self.z * v.z(),
        }
    }
}

// multiply with a f64 
impl Mul<f64> for Vec3 {
    type Output = Vec3;

    fn mul(self, q: f64) -> Vec3 {
        Vec3{
            x: self.x * q,
            y: self.y * q,
            z: self.z * q,
        }
    }
}

// divide by a Vec3
impl Div<Vec3> for Vec3 {
    type Output = Vec3;

    fn div(self, q: Vec3) -> Vec3 {
        Vec3 {
            x: self.x / q.x(),
            y: self.y / q.y(),
            z: self.z / q.z(),
        }
    }
}

// Devide by a float
impl Div<f64> for Vec3 {
    type Output = Vec3;

    fn div(self, q: f64) -> Vec3 {
        Vec3 {
            x: self.x / q,
            y: self.y / q,
            z: self.z / q,
        }
    }
}

impl PartialEq for Vec3 {
    fn eq(&self, w: &Vec3) -> bool {
        self.x == w.x() && self.y == w.y() && self.z == w.z()
    }
}

#[cfg(test)]
mod tests {
    use assert_approx_eq::assert_approx_eq;
    use crate::vec3::Vec3;

    // Test creating a new vector with three rows, taking floats as an argument
    #[test_log::test]
    fn test_create() {
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
    fn test_neg() {
        let p = Vec3::new(0.1, 0.2, 0.3);
        let q = -p;
        assert_approx_eq!(q.x(), -0.1);
        assert_approx_eq!(q.y(), -0.2);
        assert_approx_eq!(q.z(), -0.3);
    }

    #[test_log::test]
    fn test_add() {
        let p = Vec3::new(0.1, 0.2, 0.3);
        let q = Vec3::new(0.2, 0.3, 0.4);
        let r = p + q;
        assert_approx_eq!(r.x(), 0.3);
        assert_approx_eq!(r.y(), 0.5);
        assert_approx_eq!(r.z(), 0.7);
    }

    #[test_log::test]
    fn test_sub() {
        let p = Vec3::new(0.1, 0.2, 0.3);
        let q = Vec3::new(0.2, 0.3, 0.4);
        let r = p - q;
        assert_approx_eq!(r.x(), -0.1);
        assert_approx_eq!(r.y(), -0.1);
        assert_approx_eq!(r.z(), -0.1);
    }

    #[test_log::test]
    fn test_mul() {
        let p = Vec3::new(0.1, 0.2, 0.3);
        let q = Vec3::new(0.2, 0.3, 0.4);
        let r = p * q;
        assert_approx_eq!(r.x(), 0.02);
        assert_approx_eq!(r.y(), 0.06);
        assert_approx_eq!(r.z(), 0.12);
    }
    
    #[test_log::test]
    fn test_div() {
        let p = Vec3::new(0.1, 0.2, 0.3);
        let q = Vec3::new(0.2, 0.3, 0.4);
        let r = p / q;
        assert_approx_eq!(r.x(), 0.5);
        assert_approx_eq!(r.y(), 0.6666666666666666);
        assert_approx_eq!(r.z(), 0.3 / 0.4);
    }    

    // Test vector distance
    #[test_log::test]
    fn test_distance() {
        let p: Vec3 = Vec3::new(1.0, 0.0, 5.0);
        let q: Vec3 = Vec3::new(0.0, 2.0, 4.0);
        assert_approx_eq!(p.distance(&q), 2.44948974278);
        
        let w: Vec3 = Vec3::new(2.0, 3.0, 5.0);
        let v: Vec3 = Vec3::new(2.0, 0.0, 9.0);
        assert_approx_eq!(w.distance(&v), 5.0);        
    }

    // Test vector length
    #[test_log::test]
    fn test_length() {
        let p: Vec3 = Vec3::new(3.0, 1.0, 2.0);
        assert_approx_eq!(p.length(), 3.741657386773941);
        
        let q: Vec3 = Vec3::new(2.0, 3.0, 6.0);
        assert_approx_eq!(q.length(), 7.0); 
        
        let z: Vec3 = Vec3::new(0.0,0.0,0.0);
        assert_approx_eq!(z.length(), 0.0);     
    }

    #[test_log::test]
    fn test_dot() {
        let p = Vec3::new(0.1, 0.2, 0.3);
        let q = Vec3::new(0.2, 0.3, 0.4);
        assert_approx_eq!(p.dot(&q), 0.2);
    }

    // Test vector normalization
    #[test_log::test]
    fn test_normalize() {
        let q: Vec3 = Vec3::new(0.0, 2.0, 0.0).normalized();
        assert_approx_eq!(q.x(),0.0);
        assert_approx_eq!(q.y(),1.0);
        assert_approx_eq!(q.z(),0.0);
    
        let p: Vec3 = Vec3::new(3.0, 1.0, 2.0).normalized();
        assert_approx_eq!(p.x(),0.801783725737273);
        assert_approx_eq!(p.y(),0.267261241912424);
        assert_approx_eq!(p.z(),0.534522483824848);
    }
}