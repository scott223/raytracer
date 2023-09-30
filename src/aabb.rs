use crate::interval::Interval;
use crate::ray::Ray;
use crate::vec3::Vec3;
use std::mem;

// axis aligned bounding box

#[derive(Debug, Clone, Copy)]
pub struct Aabb {
    x: Interval,
    y: Interval,
    z: Interval,
}

impl Default for Aabb {
    fn default() -> Self  {
        Aabb {
            x: Interval { interval_min: 0.0, interval_max:0.0 },
            y: Interval { interval_min: 0.0, interval_max:0.0 },
            z: Interval { interval_min: 0.0, interval_max:0.0 },
        }
    }
}

impl Aabb {
    
    pub fn new_from_intervals(x: Interval, y: Interval, z: Interval) -> Self {
        Aabb { x, y, z }
    }

    pub fn new_from_points(p: Vec3, q: Vec3) -> Self {
        Aabb {
            x: Interval::new(q.x().min(p.x()), q.x().max(p.x())),
            y: Interval::new(q.y().min(p.y()), q.y().max(p.y())),
            z: Interval::new(q.z().min(p.z()), q.z().max(p.z())),
        }
    }

    pub fn new_from_aabbs(a: Aabb, b: Aabb) -> Self {
        Aabb {
            x: Interval::new_from_intervals(a.x, b.x),
            y: Interval::new_from_intervals(a.y, b.y),
            z: Interval::new_from_intervals(a.z, b.z),
        }
    }

    // return an AABB that has no side narrower than some delta, padding if necessary
    pub fn pad(&self) -> Self { 
        let delta: f64 = 0.0001;

        let new_x: Interval = if self.x.size() >= delta { self.x } else {self.x.expand(delta)};
        let new_y: Interval = if self.y.size() >= delta { self.y } else {self.y.expand(delta)};
        let new_z: Interval = if self.z.size() >= delta { self.z } else {self.z.expand(delta)};

        Aabb {
            x: new_x,
            y: new_y,
            z: new_z,
        }
    }

    pub fn axis(&self, n: usize) -> Interval {
        match n {
            0 => self.x,
            1 => self.y,
            2 => self.z,
            _ => panic!("axis out of bounds"),
        }
    }

    // checks if we have a hit with the aabb, in a given interval
    pub fn hit(&self, ray: &Ray, ray_t: &mut Interval) -> bool {
        for a in 0..3 as usize {
            let mut t0: f64 = (self.axis(a).interval_min - ray.origin.axis(a)) * ray.inv_dir.axis(a);
            let mut t1: f64 = (self.axis(a).interval_max - ray.origin.axis(a)) * ray.inv_dir.axis(a);

            // we need to swap t0 and t1
            if ray.inv_dir.axis(a) < 0.0 {
                mem::swap(&mut t0, &mut t1);
            }
            //println!("inv_d: {}, t0: {}, t1: {}", inv_d, t0, t1);

            let mut check_int: Interval = ray_t.clone();

            if t0 > check_int.interval_min {
                check_int.interval_min = t0
            };
            if t1 < check_int.interval_max {
                check_int.interval_max = t1
            };

            // no hit wit the aabb
            if check_int.interval_max <= check_int.interval_min {
                return false;
            }
        }

        // we have a hit with the aabb, so return true
        return true;
    }
}

#[cfg(test)]
mod tests {
    use crate::aabb::Aabb;
    use crate::interval::Interval;
    use crate::ray::Ray;
    use crate::vec3::Vec3;
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
}
