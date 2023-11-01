use std::f64::consts::PI;
use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Sub};

use rand::Rng;

use rand_distr::Distribution;
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct Vec3 {
    x: f64,
    y: f64,
    z: f64,
}

impl Vec3 {
    // Creating a new vector
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Vec3 { x, y, z }
    }

    // creates a random vector between a min and a max value
    pub fn new_random(val_min: f64, val_max: f64, rng: &mut impl Rng) -> Self {
        //let mut small_rng = SmallRng::from_entropy();

        Vec3 {
            x: rng.gen_range(val_min..val_max),
            y: rng.gen_range(val_min..val_max),
            z: rng.gen_range(val_min..val_max),
        }
    }

    // creates a vector within a unit sphere (lenght < 1.0)
    pub fn new_random_unit_sphere(rng: &mut impl Rng) -> Self {
        let [x, y, z] = rand_distr::UnitSphere.sample(rng);
        Self { x, y, z }
    }

    // creates a random point on a unit disk (length <1.0, z = 0)
    pub fn new_random_in_unit_disk(rng: &mut impl Rng) -> Self {
        loop {
            let v = Vec3::new(rng.gen_range(-1.0..1.0), rng.gen_range(-1.0..1.0), 0.0);
            if v.length_squared() < 1.0 {
                return v;
            }
        }
    }

    // creates a normalized random vector in the unit sphere (initial lenght < 1.0)
    pub fn new_random_unit_vector(rng: &mut impl Rng) -> Self {
        Vec3::new_random_unit_sphere(rng) //this is already normalized
    }

    //creates a normalized random vector (on unit sphere) on the same hemisphere as given vector
    pub fn new_random_on_hemisphere(normal: &Vec3, rng: &mut impl Rng) -> Self {
        let rv = Vec3::new_random_unit_vector(rng);
        if rv.dot(normal) > 0.0 {
            //in the same hemisphere
            rv
        } else {
            -rv
        }
    }

    //creates a new random cosine pdf weigthed vector
    pub fn random_cosine_direction(rng: &mut impl Rng) -> Self {
        let r1: f64 = rng.gen_range(0.0..1.0);
        let r2: f64 = rng.gen_range(0.0..1.0);

        let phi: f64 = 2. * PI * r1;
        let x: f64 = phi.cos() * r2.sqrt();
        let y: f64 = phi.sin() * r2.sqrt();
        let z: f64 = (1. - r2).sqrt();

        Vec3::new(x, y, z)
    }

    // checks if the vector is almost zero
    pub fn near_zero(&self) -> bool {
        self.x.abs() < f64::EPSILON && self.y.abs() < f64::EPSILON && self.z.abs() < f64::EPSILON
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

    pub fn axis(&self, n: usize) -> f64 {
        match n {
            1 => self.y,
            2 => self.z,
            _ => self.x,
        }
    }

    // fn distance
    // calculate the Euclidean distance between two points
    // d(p, w) = sqrt((px - wx)^2 + (py - wy)^2 + (pz - wz)^2)
    // TODO add benchmark and replace with x * x to see if its
    pub fn distance(&self, w: &Vec3) -> f64 {
        let dx = (self.x - w.x()).powf(2.0);
        let dy = f64::powf(self.y - w.y(), 2.0);
        let dz = f64::powf(self.z - w.z(), 2.0);
        (dx + dy + dz).sqrt()
    }

    // fn length
    // calculate the vector (length) between point p and the origin (0,0,0)
    pub fn length(&self) -> f64 {
        self.distance(&Vec3::new(0.0, 0.0, 0.0))
    }

    // fn lenght
    // calculates the length, squared
    pub fn length_squared(&self) -> f64 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    // fn abs
    // returns absolute value of vector, which is the same as its length (distance from 0.0)
    pub fn abs(&self) -> f64 {
        self.length()
    }

    // vector dot multiplication of self with a given vector
    pub fn dot(&self, w: &Vec3) -> f64 {
        self.x * w.x + self.y * w.y + self.z * w.z
    }

    // vector cross product
    pub fn cross(&self, other: &Vec3) -> Vec3 {
        Vec3::new(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )
    }

    // fn normalized
    // return the normalized vector (= unit vector)
    pub fn normalized(&self) -> Self {
        let l = self.length();
        Vec3::new(self.x() / l, self.y() / l, self.z() / l)
    }

    // returns a vector that is made from the smallest components of the two vectors
    pub fn min(&self, other: Vec3) -> Self {
        Vec3 {
            x: self.x().min(other.x()),
            y: self.y().min(other.y()),
            z: self.z().min(other.z()),
        }
    }

    // returns a vector that is made from the smallest components of the two vectors
    pub fn max(&self, other: Vec3) -> Self {
        Vec3 {
            x: self.x().max(other.x()),
            y: self.y().max(other.y()),
            z: self.z().max(other.z()),
        }
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
        Vec3 {
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
        Vec3 {
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

// Divide by a f64
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

// == operator (is equal)
impl PartialEq for Vec3 {
    fn eq(&self, w: &Vec3) -> bool {
        self.x == w.x() && self.y == w.y() && self.z == w.z()
    }
}

#[cfg(test)]
mod tests {
    use crate::vec3::Vec3;
    use assert_approx_eq::assert_approx_eq;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    // Test creating a new vector with three rows, taking floats as an argument
    #[test_log::test]
    fn test_new() {
        let p: Vec3 = Vec3::new(0.1, 0.2, 0.3);
        assert_approx_eq!(p.x(), 0.1);
        assert_approx_eq!(p.y(), 0.2);
        assert_approx_eq!(p.z(), 0.3);

        let q: Vec3 = Vec3::new(0.2, 0.3, 0.4);
        assert_approx_eq!(q.x(), 0.2);
        assert_approx_eq!(q.y(), 0.3);
        assert_approx_eq!(q.z(), 0.4);
    }

    #[test_log::test]
    fn test_axis() {
        let p: Vec3 = Vec3::new(0.1, 0.2, 0.3);

        assert_approx_eq!(p.axis(0), 0.1);
        assert_approx_eq!(p.axis(1), 0.2);
        assert_approx_eq!(p.axis(2), 0.3);
    }

    #[test_log::test]
    fn test_new_random() {
        let p: Vec3 = Vec3::new_random(-1.0, 1.0, &mut StdRng::seed_from_u64(222));

        assert!(p.x() < 1.01 && p.x() > -1.01);
        assert!(p.y() < 1.01 && p.y() > -1.01);
        assert!(p.z() < 1.01 && p.z() > -1.01);
    }

    #[test_log::test]
    fn test_new_random_unit_sphere() {
        let p: Vec3 = Vec3::new_random_unit_sphere(&mut StdRng::seed_from_u64(223));

        assert!(p.x() < 1.01 && p.x() > -1.01);
        assert!(p.y() < 1.01 && p.y() > -1.01);
        assert!(p.z() < 1.01 && p.z() > -1.01);
        assert!(p.length() < 1.0);
    }

    #[test_log::test]
    fn test_new_random_in_unit_disk() {
        let p: Vec3 = Vec3::new_random_in_unit_disk(&mut StdRng::seed_from_u64(223));

        assert!(p.x() < 1.01 && p.x() > -1.01);
        assert!(p.y() < 1.01 && p.y() > -1.01);
        assert!(p.z() == 0.0);
        assert!(p.length() < 1.0);
    }

    #[test_log::test]
    fn test_new_random_unit_vector() {
        let p: Vec3 = Vec3::new_random_unit_vector(&mut StdRng::seed_from_u64(223));

        assert!(p.x() < 1.01 && p.x() > -1.01);
        assert!(p.y() < 1.01 && p.y() > -1.01);
        assert!(p.z() < 1.01 && p.z() > -1.01);
        assert_approx_eq!(p.length(), 1.0);
    }

    #[test_log::test]
    fn test_new_random_hemisphere() {
        let n: Vec3 = Vec3::new(0.0, 1.0, 0.0);
        let p: Vec3 = Vec3::new_random_on_hemisphere(&n, &mut StdRng::seed_from_u64(222)); // this is a unit vector with lenght = 1

        assert!(p.x() < 1.01 && p.x() > -1.01);
        assert!(p.y() < 1.01 && p.y() > -1.01);
        assert!(p.z() < 1.01 && p.z() > -1.01);

        assert!(p.dot(&n) > 0.0);
    }

    #[test_log::test]
    fn test_min() {
        let n: Vec3 = Vec3::new(1.0, 2.0, 3.0);
        let p: Vec3 = Vec3::new(4.0, 3.0, 2.0);

        assert_approx_eq!(n.min(p), Vec3::new(1.0, 2.0, 2.0));
        //assert!(p.near_zero());
    }

    #[test_log::test]
    fn test_max() {
        let n: Vec3 = Vec3::new(1.0, 2.0, 3.0);
        let p: Vec3 = Vec3::new(4.0, 3.0, 2.0);

        assert_approx_eq!(n.max(p), Vec3::new(4.0, 3.0, 3.0));
        //assert!(p.near_zero());
    }

    #[test_log::test]
    fn test_near_zero() {
        let n: Vec3 = Vec3::new(0.0, 0.001, 0.0);
        let p: Vec3 = Vec3::new(0.0, 0.0, 0.0);

        assert!(!n.near_zero());
        assert!(p.near_zero());
    }

    #[test_log::test]
    fn test_neg() {
        let p: Vec3 = Vec3::new(0.1, 0.2, 0.3);
        let q: Vec3 = -p;
        assert_approx_eq!(q.x(), -0.1);
        assert_approx_eq!(q.y(), -0.2);
        assert_approx_eq!(q.z(), -0.3);
    }

    #[test_log::test]
    fn test_add() {
        let p: Vec3 = Vec3::new(0.1, 0.2, 0.3);
        let q: Vec3 = Vec3::new(0.2, 0.3, 0.4);
        let r: Vec3 = p + q;
        assert_approx_eq!(r.x(), 0.3);
        assert_approx_eq!(r.y(), 0.5);
        assert_approx_eq!(r.z(), 0.7);
    }

    #[test_log::test]
    fn test_sub() {
        let p: Vec3 = Vec3::new(0.1, 0.2, 0.3);
        let q: Vec3 = Vec3::new(0.2, 0.3, 0.4);
        let r: Vec3 = p - q;
        assert_approx_eq!(r.x(), -0.1);
        assert_approx_eq!(r.y(), -0.1);
        assert_approx_eq!(r.z(), -0.1);
    }

    #[test_log::test]
    fn test_mul() {
        let p: Vec3 = Vec3::new(0.1, 0.2, 0.3);
        let q: Vec3 = Vec3::new(0.2, 0.3, 0.4);
        let r: Vec3 = p * q;
        assert_approx_eq!(r.x(), 0.02);
        assert_approx_eq!(r.y(), 0.06);
        assert_approx_eq!(r.z(), 0.12);
    }

    #[test_log::test]
    fn test_div() {
        let p: Vec3 = Vec3::new(0.1, 0.2, 0.3);
        let q: Vec3 = Vec3::new(0.2, 0.3, 0.4);
        let r: Vec3 = p / q;
        assert_approx_eq!(r.x(), 0.5);
        assert_approx_eq!(r.y(), 0.6666666666666666);
        assert_approx_eq!(r.z(), 0.3 / 0.4);
    }

    #[test_log::test]
    fn test_distance() {
        let p: Vec3 = Vec3::new(1.0, 0.0, 5.0);
        let q: Vec3 = Vec3::new(0.0, 2.0, 4.0);
        assert_approx_eq!(p.distance(&q), 2.44948974278);

        let w: Vec3 = Vec3::new(2.0, 3.0, 5.0);
        let v: Vec3 = Vec3::new(2.0, 0.0, 9.0);
        assert_approx_eq!(w.distance(&v), 5.0);
    }

    #[test_log::test]
    fn test_length() {
        let p: Vec3 = Vec3::new(3.0, 1.0, 2.0);
        assert_approx_eq!(p.length(), 3.741657386773941);

        let q: Vec3 = Vec3::new(2.0, 3.0, 6.0);
        assert_approx_eq!(q.length(), 7.0);

        let z: Vec3 = Vec3::new(0.0, 0.0, 0.0);
        assert_approx_eq!(z.length(), 0.0);
    }

    #[test_log::test]
    fn test_length_squared() {
        let p: Vec3 = Vec3::new(3.0, 1.0, 2.0);
        assert_approx_eq!(p.length_squared(), 14.0);

        let q: Vec3 = Vec3::new(2.0, 3.0, 6.0);
        assert_approx_eq!(q.length_squared(), 49.0);

        let z: Vec3 = Vec3::new(0.0, 0.0, 0.0);
        assert_approx_eq!(z.length_squared(), 0.0);
    }

    #[test_log::test]
    fn test_dot() {
        let p: Vec3 = Vec3::new(0.1, 0.2, 0.3);
        let q: Vec3 = Vec3::new(0.2, 0.3, 0.4);
        assert_approx_eq!(p.dot(&q), 0.2);
    }

    #[test_log::test]
    fn test_cross() {
        let p: Vec3 = Vec3::new(1.0, 0.0, 0.0);
        let q: Vec3 = Vec3::new(3.0, 2.0, 4.0);
        assert_approx_eq!(p.cross(&q), Vec3::new(0.0, -4.0, 2.0));
    }

    // Test vector normalization
    #[test_log::test]
    fn test_normalize() {
        let q: Vec3 = Vec3::new(0.0, 2.0, 0.0).normalized();
        assert_approx_eq!(q.x(), 0.0);
        assert_approx_eq!(q.y(), 1.0);
        assert_approx_eq!(q.z(), 0.0);

        let p: Vec3 = Vec3::new(3.0, 1.0, 2.0).normalized();
        assert_approx_eq!(p.x(), 0.801783725737273);
        assert_approx_eq!(p.y(), 0.267261241912424);
        assert_approx_eq!(p.z(), 0.534522483824848);
    }
}
