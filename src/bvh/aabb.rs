use crate::elements::Element;
use crate::elements::Hittable;
use crate::linalg::Vec3;
use crate::render::Axis;
use crate::render::Interval;
use crate::render::Ray;

use std::fmt;
use std::ops::Index;

/// Axis Aligned Bounding Box struct
#[derive(Debug, Clone, Copy)]
pub struct Aabb {
    x: Interval,
    y: Interval,
    z: Interval,
    pub min: Vec3,
    pub max: Vec3,
    pub centroid: Vec3,
}

impl Default for Aabb {
    fn default() -> Self {
        Aabb {
            x: Interval {
                interval_min: 0.0,
                interval_max: 0.0,
            },
            y: Interval {
                interval_min: 0.0,
                interval_max: 0.0,
            },
            z: Interval {
                interval_min: 0.0,
                interval_max: 0.0,
            },
            min: Vec3::new(0.0, 0.0, 0.0),
            max: Vec3::new(0.0, 0.0, 0.0),
            centroid: Vec3::new(0.0, 0.0, 0.0),
        }
    }
}

impl Aabb {
    pub fn new_from_intervals(x: Interval, y: Interval, z: Interval) -> Self {
        let mut bbox = Aabb {
            x,
            y,
            z,
            min: Vec3::new(x.interval_min, y.interval_min, z.interval_min),
            max: Vec3::new(x.interval_max, y.interval_max, z.interval_max),
            centroid: Vec3::new(0.0, 0.0, 0.0),
        };

        bbox.centroid = calculate_centroid(bbox);
        bbox
    }

    pub fn new_from_points(p: Vec3, q: Vec3) -> Self {
        let x: Interval = Interval::new(q.x().min(p.x()), q.x().max(p.x()));
        let y: Interval = Interval::new(q.y().min(p.y()), q.y().max(p.y()));
        let z: Interval = Interval::new(q.z().min(p.z()), q.z().max(p.z()));

        let mut bbox = Aabb {
            x,
            y,
            z,
            min: Vec3::new(x.interval_min, y.interval_min, z.interval_min),
            max: Vec3::new(x.interval_max, y.interval_max, z.interval_max),
            centroid: Vec3::new(0.0, 0.0, 0.0),
        };

        bbox.centroid = calculate_centroid(bbox);
        bbox
    }

    pub fn new_from_point(p: Vec3) -> Self {
        Aabb::new_from_points(p, p)
    }

    pub fn new_from_aabbs(a: Aabb, b: Aabb) -> Self {
        let x: Interval = Interval::new_from_intervals(a.x, b.x);
        let y: Interval = Interval::new_from_intervals(a.y, b.y);
        let z: Interval = Interval::new_from_intervals(a.z, b.z);

        let mut bbox = Aabb {
            x,
            y,
            z,
            min: Vec3::new(x.interval_min, y.interval_min, z.interval_min),
            max: Vec3::new(x.interval_max, y.interval_max, z.interval_max),
            centroid: Vec3::new(0.0, 0.0, 0.0),
        };

        bbox.centroid = calculate_centroid(bbox);
        bbox
    }

    pub fn new_from_objects(objects: &[Element], indices: &Vec<usize>) -> Self {
        // start met bbox from first element

        let mut bbox = objects[indices[0]].bounding_box();

        if indices.len() > 1 {
            // walk through each object and add the bounding box
            for index in indices {
                bbox = Aabb::new_from_aabbs(bbox, objects[*index].bounding_box());
            }
        }

        bbox
    }

    /// Returns the surface area of this Aabb
    pub fn area(&self) -> f64 {
        (self.max.x() - self.min.x()) * (self.max.y() - self.min.y())
            + (self.max.y() - self.min.y()) * (self.max.z() - self.min.z())
            + (self.max.z() - self.min.z()) * (self.max.x() - self.min.x())
    }

    /// Returns an AABB that has no side narrower than some delta, padding if necessary.
    pub fn pad(&self) -> Self {
        let delta: f64 = 0.0001;

        let new_x: Interval = if self.x.size() >= delta {
            self.x
        } else {
            self.x.expand(delta)
        };
        let new_y: Interval = if self.y.size() >= delta {
            self.y
        } else {
            self.y.expand(delta)
        };
        let new_z: Interval = if self.z.size() >= delta {
            self.z
        } else {
            self.z.expand(delta)
        };

        let mut bbox = Aabb {
            x: new_x,
            y: new_y,
            z: new_z,
            min: Vec3::new(new_x.interval_min, new_y.interval_min, new_z.interval_min),
            max: Vec3::new(new_x.interval_max, new_y.interval_max, new_z.interval_max),
            centroid: Vec3::new(0.0, 0.0, 0.0),
        };

        bbox.centroid = calculate_centroid(bbox);
        bbox
    }

    #[inline(always)]
    pub fn axis(&self, n: usize) -> Interval {
        match n {
            0 => self.x,
            1 => self.y,
            2 => self.z,
            _ => panic!("axis out of bounds"),
        }
    }

    #[inline(always)]
    pub fn largest_axis(&self) -> Axis {
        let size = self.max - self.min;

        if size.x() > size.y() && size.x() > size.z() {
            Axis::X
        } else if size.y() > size.z() {
            Axis::Y
        } else {
            Axis::Z
        }
    }

    /// Returns the position of the centriod of the Aabb.
    #[inline(always)]
    pub fn centroid(&self) -> Vec3 {
        calculate_centroid(*self)
    }

    /// Returns if we have a hit with the aabb, in a given interval
    ///
    /// Source: <https://docs.rs/bvh/latest/src/bvh/ray.rs.html#168-188>

    pub fn hit(&self, ray: &Ray, _ray_t: &mut Interval) -> bool {
        let mut ray_min = (self[ray.sign_x].x() - ray.origin.x()) * ray.inv_direction.x();
        let mut ray_max = (self[1 - ray.sign_x].x() - ray.origin.x()) * ray.inv_direction.x();

        let y_min = (self[ray.sign_y].y() - ray.origin.y()) * ray.inv_direction.y();
        let y_max = (self[1 - ray.sign_y].y() - ray.origin.y()) * ray.inv_direction.y();

        ray_min = ray_min.max(y_min);
        ray_max = ray_max.min(y_max);

        let z_min = (self[ray.sign_z].z() - ray.origin.z()) * ray.inv_direction.z();
        let z_max = (self[1 - ray.sign_z].z() - ray.origin.z()) * ray.inv_direction.z();

        ray_min = ray_min.max(z_min);
        ray_max = ray_max.min(z_max);

        ray_min.max(0.0) <= ray_max
    }

    #[inline(always)]
    pub fn grow<T: Hittable>(mut self, b: &T) {
        self.x = Interval::new_from_intervals(self.x, b.bounding_box().x);
        self.y = Interval::new_from_intervals(self.y, b.bounding_box().y);
        self.z = Interval::new_from_intervals(self.z, b.bounding_box().z);
    }
}

pub fn calculate_centroid(bbox: Aabb) -> Vec3 {
    Vec3::new(
        bbox.min.x() + (bbox.max.x() - bbox.min.x()) / 2.,
        bbox.min.y() + (bbox.max.y() - bbox.min.y()) / 2.,
        bbox.min.z() + (bbox.max.z() - bbox.min.z()) / 2.,
    )
}

/// Display trait for Aabb
impl fmt::Display for Aabb {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Min bound: {}; Max bound: {}", self.min, self.max)
    }
}

/// Make [`Aabb`]s indexable. `aabb[0]` gives a reference to the minimum bound.
/// All other indices return a reference to the maximum bound.
impl Index<usize> for Aabb {
    type Output = Vec3;

    fn index(&self, index: usize) -> &Vec3 {
        if index == 0 {
            &self.min
        } else {
            &self.max
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::bvh::aabb::Aabb;
    use crate::linalg::Vec3;
    use crate::render::Axis;
    use crate::render::Interval;
    use crate::render::Ray;
    use assert_approx_eq::assert_approx_eq;

    // Test creating a new aabb, taking intervals as an argument
    #[test_log::test]
    fn test_new_from_interval() {
        let x: Interval = Interval::new(0.1, 0.2);
        let y: Interval = Interval::new(0.3, 0.4);
        let z: Interval = Interval::new(0.2, 0.3);

        let aabb_p = Aabb::new_from_intervals(x, y, z);

        assert_approx_eq!(aabb_p.x.interval_min, x.interval_min);
    }

    // Test creating a new aabb, taking points as an argument
    #[test_log::test]
    fn test_new_from_points() {
        let p: Vec3 = Vec3::new(0.1, 0.1, 0.1);
        let q: Vec3 = Vec3::new(0.2, 0.0, 0.2);

        let aabb_p: Aabb = Aabb::new_from_points(p, q);

        assert_approx_eq!(aabb_p.x.interval_min, 0.1);
        assert_approx_eq!(aabb_p.y.interval_min, 0.0);
        assert_approx_eq!(aabb_p.z.interval_min, 0.1);

        assert_approx_eq!(aabb_p.x.interval_max, 0.2);
        assert_approx_eq!(aabb_p.y.interval_max, 0.1);
        assert_approx_eq!(aabb_p.z.interval_max, 0.2);
    }

    #[test_log::test]
    fn test_new_from_aabbs() {
        let p: Vec3 = Vec3::new(1.0, 1.0, 1.0);
        let q: Vec3 = Vec3::new(2.0, 2.0, 2.0);

        let a: Aabb = Aabb::new_from_points(p, q);

        let w: Vec3 = Vec3::new(3.0, 3.0, 3.0);
        let v: Vec3 = Vec3::new(4.0, 4.0, 4.0);

        let b: Aabb = Aabb::new_from_points(w, v);

        let aabb_box: Aabb = Aabb::new_from_aabbs(a, b);

        assert_approx_eq!(aabb_box.x.interval_min, 1.0);
        assert_approx_eq!(aabb_box.y.interval_min, 1.0);
        assert_approx_eq!(aabb_box.z.interval_min, 1.0);

        assert_approx_eq!(aabb_box.x.interval_max, 4.0);
        assert_approx_eq!(aabb_box.y.interval_max, 4.0);
        assert_approx_eq!(aabb_box.z.interval_max, 4.0);
    }

    #[test_log::test]
    fn test_axis() {
        let p: Vec3 = Vec3::new(0.1, 0.1, 0.1);
        let q: Vec3 = Vec3::new(0.2, 0.0, 0.2);

        let aabb_p: Aabb = Aabb::new_from_points(p, q);

        assert_approx_eq!(aabb_p.axis(0).interval_min, 0.1);
        assert_approx_eq!(aabb_p.axis(0).interval_max, 0.2);

        assert_approx_eq!(aabb_p.axis(1).interval_min, 0.0);
        assert_approx_eq!(aabb_p.axis(1).interval_max, 0.1);

        assert_approx_eq!(aabb_p.axis(2).interval_min, 0.1);
        assert_approx_eq!(aabb_p.axis(2).interval_max, 0.2);
    }

    #[test_log::test]
    fn test_hit() {
        let p: Vec3 = Vec3::new(1.0, 1.0, 1.0);
        let q: Vec3 = Vec3::new(2.0, 2.0, 2.0);

        let aabb_p: Aabb = Aabb::new_from_points(p, q);
        let ray: Ray = Ray::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 1.0, 1.0));
        let mut int: Interval = Interval::new(0.0, 20.0);

        assert_eq!(aabb_p.hit(&ray, &mut int), true);

        let ray_two: Ray = Ray::new(Vec3::new(5.0, 5.0, 5.0), Vec3::new(6.0, 6.0, 6.0));

        assert_eq!(aabb_p.hit(&ray_two, &mut int), false);
    }

    #[test_log::test]
    fn test_centroid() {
        let p: Vec3 = Vec3::new(1.0, 1.0, 1.0);
        let q: Vec3 = Vec3::new(2.0, 2.0, 2.0);

        let aabb_p: Aabb = Aabb::new_from_points(p, q);

        let centroid = aabb_p.centroid();

        assert_eq!(centroid.x(), 1.5);
        assert_eq!(centroid.y(), 1.5);
        assert_eq!(centroid.z(), 1.5);
    }

    #[test_log::test]
    fn test_longest_axis() {
        let p: Vec3 = Vec3::new(1.0, 1.0, 1.0);
        let q: Vec3 = Vec3::new(2.0, 3.0, 2.0);

        let aabb_p: Aabb = Aabb::new_from_points(p, q);

        let largest_axis = aabb_p.largest_axis();

        assert_eq!(largest_axis, Axis::Y);

        let p: Vec3 = Vec3::new(1.0, 1.0, 1.0);
        let q: Vec3 = Vec3::new(2.0, 3.0, 5.0);

        let aabb_p: Aabb = Aabb::new_from_points(p, q);

        let largest_axis = aabb_p.largest_axis();

        assert_eq!(largest_axis, Axis::Z);
    }

    #[test_log::test]
    fn test_surface_area() {
        let p: Vec3 = Vec3::new(1.0, 1.0, 1.0);
        let q: Vec3 = Vec3::new(2.0, 3.0, 2.0);

        let aabb: Aabb = Aabb::new_from_points(p, q);

        let area = aabb.area();

        assert_eq!(area, (1.*2.)+(2.*1.)+(1.*1.));

        let p: Vec3 = Vec3::new(-2.0, 1.0, -2.0);
        let q: Vec3 = Vec3::new(2.0, 5.0, 2.0);

        let aabb: Aabb = Aabb::new_from_points(p, q);

        let area = aabb.area();

        assert_eq!(area, (4.*4.)+(4.*4.)+(4.*4.));
    }
}
