use crate::element::HitRecord;
use crate::element::Hittable;
use crate::ray::Ray;
use crate::vec3::Vec3;

#[derive(Debug, Clone, Copy)]
pub struct Plane {
    pub origin: Vec3,
    pub normal: Vec3,
    pub color: Vec3,
}

impl Plane {
    // Creating a new Plane with an origin and a normal (planes have infitine size)
    pub fn new(origin: Vec3, normal: Vec3) -> Self {
        Plane {
            origin: origin,
            normal: normal.normalized(),
            color: Vec3::new(30.0,30.0,30.0),
        }
    }
}

impl Hittable for Plane {
    // finding the hits for a given ray
    // https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-plane-and-ray-disk-intersection.html

    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {

        let denom = self.normal.dot(&ray.direction);

        //log::info!("rd: {}", ray.direction);

        if denom > 1e-6 {
            let v = self.origin - ray.origin;
            let distance = v.dot(&self.normal) / denom;

            //log::info!("distance: {}", distance);

            if distance >= 0.0 {

                if distance > t_min && distance < t_max {

                    //log::info!("distance: {}", distance);

                    let p = ray.at(distance);
                    let hit = HitRecord {
                        t: distance,
                        normal: -self.normal, //we need a minus here to get the defraction working, not sure why.....
                        point: p,
                        color: self.color,
                    };
                    //log::info!("hit: {:?}", hit);
                    return Some(hit);
                }
            }
        }


        None // no hits found
    }
}

#[cfg(test)]
mod tests {
    use crate::plane::Plane;
    use crate::vec3::Vec3;
    use crate::ray::{Ray};
    use crate::element::{Hittable};
    use assert_approx_eq::assert_approx_eq;

    #[test_log::test]
    fn test_create_plane() {
        let p = Plane::new(Vec3::new(0.0, -2.5, 0.0), Vec3::new(0.0, -2.0, 0.0));
        assert_approx_eq!(p.origin, Vec3::new(0.0, -2.5, 0.0));
        assert_approx_eq!(p.normal, Vec3::new(0.0, -1.0, 0.0));
    }

    #[test_log::test]
    fn test_hit_plane() {
        let r = Ray::new(Vec3::new(0.0,0.0,0.0),Vec3::new(0.0, 0.0, -1.0));
        let p = Plane::new(Vec3::new(0.0, -2.5, 0.0), Vec3::new(0.0, -2.0, 0.0));

        if let Some(hit) = p.hit(&r, 0.0, f64::MAX) {
      //      assert_eq!(hit.t,2.0);
            assert_eq!(hit.color,p.color);
      //      assert_eq!(hit.point, Vec3::new(0.0,0.0,-2.0));
        } 
    }
}