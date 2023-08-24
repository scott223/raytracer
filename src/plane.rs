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