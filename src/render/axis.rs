use std::{fmt::{Display, Formatter, self}, ops::Index};

use crate::linalg::Vec3;


/// An `Axis` in a three-dimensional coordinate system.
/// Used to access `Vec3` structs via index.
///
/// # Examples
/// ```
/// use raytracer::render::Axis;
/// use raytracer::linalg::Vec3;
///
/// let mut vec = Vec3::new(1.0, 0.5, 42.0);
///
/// assert_eq!(vec[Axis::Y], 0.5);
/// ```
///
/// [`Vec3`] is also indexable using `Axis`.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum Axis {
    /// Index of the X axis.
    X = 0,

    /// Index of the Y axis.
    Y = 1,

    /// Index of the Z axis.
    Z = 2,
}

/// Display implementation for `Axis`.
impl Display for Axis {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(
            f,
            "{}",
            match *self {
                Axis::X => "x",
                Axis::Y => "y",
                Axis::Z => "z",
            }
        )
    }
}

/// Make `Vector3` indexable by `Axis`.
impl Index<Axis> for Vec3 {
    type Output = f64;

    fn index(&self, axis: Axis) -> &f64 {
        match axis {
            Axis::X => &self.x,
            Axis::Y => &self.y,
            Axis::Z => &self.z,
        }
    }
}