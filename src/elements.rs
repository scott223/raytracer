use std::fmt::Debug;
use crate::{materials::*, aabb::Aabb};
use crate::ray::Ray;
use crate::vec3::Vec3;
use crate::interval::Interval;

use serde::{Deserialize, Serialize};

// hitrecord gets returned on a hit, containting the point on the ray, the point in the global coordinate system, the normal and the material for the hit
#[derive(Debug)]
pub struct HitRecord {
    pub t: f64,
    pub point: Vec3,
    pub normal: Vec3,
    pub front_face: bool,
    pub material: Material,
}

// hittable trait defines the hit function for each element
pub trait Hittable {
    fn hit(&self, ray: &Ray, ray_t: &mut Interval) -> Option<HitRecord>;
    fn bounding_box(&self) -> Aabb;
}

// enum for all the different elements, simplified JSON representation
#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub enum JSONElement {
    JSONSphere(JSONSphere),
    JSONQuad(JSONQuad),
    //Triangle(Triangle),
}


// enum for all the different elements
#[derive(Debug, Clone, Copy)]
pub enum Element {
    Sphere(Sphere),
    Quad(Quad),
    Triangle(Triangle),
}

// matching the hit function with the hittable trait for each type of element
impl Hittable for Element {
    fn hit(&self, ray: &Ray, ray_t: &mut Interval) -> Option<HitRecord> {
        match *self {
            Element::Sphere(ref s) => s.hit(ray, ray_t),
            Element::Quad(ref q) => q.hit(ray, ray_t),
            Element::Triangle(ref t) => t.hit(ray, ray_t),
 //           Element::Plane(ref p) => p.hit(ray, t_min, t_max),
        }
    }

    fn bounding_box(&self) -> Aabb {
        match *self {
            Element::Sphere(ref s) => s.bounding_box(),
            Element::Quad(ref q) => q.bounding_box(),
            Element::Triangle(ref t) => t.bounding_box(),
        }
    }
}

//Triangle element
#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct Triangle {
    pub v0: Vec3,
    pub v1: Vec3,
    pub v2: Vec3,
    pub material: Material,
    // pub bbox: Aabb,
    // TODO initialize the bounding box at object creation!
}

impl Triangle {
    pub fn new(v0: Vec3, v1: Vec3, v2: Vec3, material: Material) -> Self {
        Triangle{
            v0,
            v1,
            v2,
            material,
            //aabb TODO define the bounding box at creation time
        }
    }
}

// https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/ray-triangle-intersection-geometric-solution.html
impl Hittable for Triangle {
    fn hit(&self, ray: &Ray, ray_t: &mut Interval) -> Option<HitRecord> {
        //compute the plane's normal
        let v0v1: Vec3 = self.v1 - self.v0;
        let v0v2: Vec3 = self.v2 - self.v0;
        let normal = v0v1.cross(&v0v2);
        // let area: f64 = normal.length();

        //step 1 - finding P
        let ndot_ray_direction = normal.dot(&ray.direction);

        // check if the ray and the plane are parallel. if they are, no hit
        if ndot_ray_direction < 0.000001 {
            return None;
        }

        // compute d parameter
        let d: f64 = (normal * -1.0).dot(&self.v0);

        // compute t
        let t = -(normal.dot(&ray.origin) + d) / ndot_ray_direction;

        // check if t is inside the interval (before the camera, and closer than an earlier hit)
        if !ray_t.contains(t) {
            return None;
        }

        // compute the intersection point
        let point: Vec3 = ray.at(t);

        // inside-outside test

        // edge 0
        let edge0: Vec3 = self.v1 - self.v0;
        let vp0: Vec3 = point - self.v0;
        let c: Vec3 = edge0.cross(&vp0);
        if normal.dot(&c) < 0.0 {
            return None; // point is on the right side
        }

        // edge 1
        let edge1: Vec3 = self.v2 - self.v1;
        let vp1: Vec3 = point - self.v1;
        let c: Vec3 = edge1.cross(&vp1);
        if normal.dot(&c) < 0.0 {
            return None; // point is on the right side
        }

        // edge 2
        let edge2: Vec3 = self.v0 - self.v2;
        let vp2: Vec3 = point - self.v2;
        let c: Vec3 = edge2.cross(&vp2);
        if normal.dot(&c) < 0.0 {
            return None; // point is on the right side
        }

        // we have a hit, so we return a hit record
        return Some(HitRecord {
            t,
            normal, //TODO check normal (should be outward normal)
            point,
            front_face: ray.direction.dot(&normal) < 0.0,
            material: self.material,
        });

    }

    // find the minimum xyz and max xyz coordinates and use for bounding box. add some padding as triangles are usually flat..
    // https://stackoverflow.com/questions/39974191/triangle-bounding-box
    fn bounding_box(&self) -> Aabb {
        Aabb::new_from_points(
            self.v0.min(self.v1.min(self.v2)), 
            self.v0.max(self.v1.max(self.v2))
        ).pad()
    }
}

// Plane element
#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct Plane {
    pub origin: Vec3,
    pub normal: Vec3,
    pub material: Material,
}

// Creating a new Plane with an origin and a normal (unit vector) (planes have infinite size)
impl Plane {
    #[allow(dead_code)]
    pub fn new(origin: Vec3, normal: Vec3, material: Material) -> Self {
        Plane {
            origin,
            normal: normal.normalized(),
            material,
        }
    }
}

// finding the hits for a given ray
// return a HitRecord, that contains the position on the ray, the point of the hit & the material
// https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-plane-and-ray-disk-intersection.html
impl Hittable for Plane {
    fn hit(&self, ray: &Ray, ray_t: &mut Interval) -> Option<HitRecord> {
        // get the denominator
        let denom = self.normal.dot(&ray.direction);

        // check if larger than zero, and positive
        if denom > 1e-6 {
            let v = self.origin - ray.origin;
            let distance = v.dot(&self.normal) / denom;

            if distance >= 0.0 && distance > ray_t.interval_min && distance < ray_t.interval_max {
                let p = ray.at(distance);
                let hit = HitRecord {
                    t: distance,
                    normal: -self.normal, //we need a minus here to get the defraction working, not sure why.....
                    point: p,
                    front_face: true,
                    material: self.material,
                };
                return Some(hit);
            }
        }
        None // no hits found
    }

    fn bounding_box(&self) -> Aabb {
        // TODO, will need to figure out how to add a bounding box for an infinite plane
        todo!();
    }
}

//Quad element

// simplified JSON version of the quad
#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct JSONQuad {
    pub Q: Vec3,
    pub u: Vec3,
    pub v: Vec3,
    pub material: Material,
}

// actual quad element
#[derive(Debug, Clone, Copy)]
pub struct Quad {
    pub Q: Vec3,
    pub u: Vec3,
    pub v: Vec3,
    pub material: Material,

    //following fields are used for every hit calculation, so we pre-calculate in the constructor
    pub n: Vec3,
    pub normal: Vec3,
    pub D: f64,
    pub w: Vec3,
    pub bbox: Aabb,
}


impl Quad {
    pub fn new_from_json_object(json_quad: JSONQuad) -> Self {
        Quad::new(json_quad.Q, json_quad.u, json_quad.v, json_quad.material)
    }
    
    // Creating a new Quad with lower left point Q and vectors u and v
    pub fn new(Q: Vec3, u: Vec3, v:Vec3, material: Material) -> Self {
        
        let n: Vec3 = u.cross(&v);
        let normal: Vec3 = n.normalized();
        let D: f64 = normal.dot(&Q);
        let w: Vec3 = n / n.dot(&n);
        
        Quad {
            Q,
            u,
            v,
            n,
            normal,
            D,
            w,
            material,
            bbox: Aabb::new_from_points(Q, Q + u + v).pad(),
        }
    }

    //given the hit point in plane coordinates, return none if it is outside the primitive
    pub fn is_interior(a: f64, b: f64) -> Option<Vec<f64>> {
        if (a < 0.0) || (1.0 < a) || (b < 0.0) || (1.0 < b) {
            return None
        }

        // we have this work around as we dont determine the u and v in the hitrecord yet
       Some(vec![a, b])
    }
}

impl Hittable for Quad {
    fn hit(&self, ray: &Ray, ray_t: &mut Interval) -> Option<HitRecord> {

        let denom = self.normal.dot(&ray.direction);

        // no hits, as the ray is parallel to the plane
        if denom.abs() < 0.000000001 {
            return None
        }

        // return false if the hit point paramater t is outside the ray interval
        let t: f64 = (self.D - self.normal.dot(&ray.origin)) / denom;
        if !ray_t.contains(t) {
            return None;
        }

        let intersection = ray.at(t);
        let planar_hitpoint_vector: Vec3 = intersection - self.Q;
        let alpha: f64 = self.w.dot(&planar_hitpoint_vector.cross(&self.v));
        let beta: f64 = self.w.dot(&self.u.cross(&planar_hitpoint_vector));

        match Quad::is_interior(alpha, beta) {
            Some(_result) => {
                // nothing needed here i think
            },
            _ => { 
                //is not in the interior, so break
                return None 
            }
        }

        // we have a hit, so we return a hit record
        return Some(HitRecord {
            t,
            normal: if ray.direction.dot(&self.normal) < 0.0 { self.normal } else { -self.normal }, //TODO check normal (should be outward normal)
            point: intersection,
            front_face: ray.direction.dot(&self.normal) < 0.0,
            material: self.material,
        });

    }
    
    fn bounding_box(&self) -> Aabb {
        //Aabb::new_from_points(self.Q, self.Q + self.u + self.v).pad()
        
        //we now calculate the bbox on initialization, so we can just return the field
        self.bbox
    }
}

// Sphere element

// simplified JSON representation of the Sphere
#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct JSONSphere {
    pub center: Vec3,
    pub radius: f64,
    pub material: Material,
}

// the actual sphere element
#[derive(Debug, Clone, Copy)]
pub struct Sphere {
    pub center: Vec3,
    pub radius: f64,
    pub material: Material,
    pub bbox: Aabb,
}

// Creating a new Sphere with a center and a radius
impl Sphere {
    
    // create a new sphere based on the JSON object
    pub fn new_from_json_object(json_sphere: JSONSphere) -> Self {
        Sphere::new(json_sphere.center, json_sphere.radius, json_sphere.material)
    }
    
    // create a new sphere
    pub fn new(center: Vec3, radius: f64, material: Material) -> Self {
        let rvec = Vec3::new(radius, radius, radius);
        let bbox = Aabb::new_from_points(center - rvec, center + rvec);
        
        Sphere {
            center,
            radius,
            material,
            bbox,
        }
    }
}

// implementing hte hittable traits for the Sphere
impl Hittable for Sphere {

    // finding the hits for a given ray
    // based on: https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-sphere-intersection.html
    fn hit(&self, ray: &Ray, ray_t: &mut Interval) -> Option<HitRecord> {
        let l = ray.origin - self.center;

        // solve the quadratic equation
        let a = ray.direction.dot(&ray.direction);
        let b = 2.0 * ray.direction.dot(&l);
        let c = l.dot(&l) - (self.radius * self.radius);

        let discr: f64 = b * b - 4.0 * a * c;

        if discr >= 0.0 {
            // there is one or two solutions

            let q: f64 = if b > 0.0 {
                -0.5 * (b + discr.sqrt())
            } else {
                -0.5 * (b - discr.sqrt())
            };
            let first_root: f64 = q / a;
            let second_root: f64 = c / q;

            if first_root < 0.0 && second_root < 0.0 {
                // only negative solutions
                return None;
            }

            let nearest_root = if first_root < second_root {
                first_root
            } else {
                second_root
            };

            if nearest_root < ray_t.interval_max && nearest_root > ray_t.interval_min {
                let p = ray.at(nearest_root); //we dont need to find the solution with the discrimant, but can just ask the ray where it was at a given t
                let outward_normal = ((p - self.center) / self.radius).normalized();
                let front_face = ray.direction.dot(&outward_normal) < 0.0;

                return Some(HitRecord {
                    t: nearest_root,
                    normal: if front_face { outward_normal } else { -outward_normal },
                    point: p,
                    front_face: front_face,
                    material: self.material,
                });
            }
        }
        None // no hits found
    }

    // construct an axis aligned bounding box Aabb for a sphere
    fn bounding_box(&self) -> Aabb {
        //let rvec = Vec3::new(self.radius, self.radius, self.radius);
        //Aabb::new_from_points(self.center - rvec, self.center + rvec)

        self.bbox
    }
}

#[cfg(test)]
mod tests {
    use crate::color::Color;
    use crate::elements::*;
    use crate::ray::Ray;
    use crate::vec3::Vec3;
    use assert_approx_eq::assert_approx_eq;

    #[test_log::test]
    fn test_create_plane() {
        let m1: Material = Material::Lambertian(Lambertian::new(Color::new(1.0, 1.0, 1.0)));
        let p: Plane = Plane::new(Vec3::new(0.0, -2.5, 0.0), Vec3::new(0.0, -2.0, 0.0), m1);
        assert_approx_eq!(p.origin, Vec3::new(0.0, -2.5, 0.0));
        assert_approx_eq!(p.normal, Vec3::new(0.0, -1.0, 0.0));
    }

    #[test_log::test]
    fn test_hit_plane() {
        let m1: Material = Material::Lambertian(Lambertian::new(Color::new(1.0, 1.0, 1.0)));
        let r: Ray = Ray::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0));
        let p: Plane = Plane::new(Vec3::new(0.0, 0.0, 5.0), Vec3::new(0.0, 0.0, 1.0), m1);

        if let Some(hit) = p.hit(&r, &mut Interval::new(0.0, f64::MAX)) {
            assert_eq!(hit.t, 5.0);
            assert_eq!(hit.point, Vec3::new(0.0, 0.0, 5.0));
        }

        let r: Ray = Ray::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, -1.0, 0.0));
        let p: Plane = Plane::new(Vec3::new(0.0, -2.5, 0.0), Vec3::new(0.0, -1.0, 0.0), m1);

        if let Some(hit) = p.hit(&r, &mut Interval::new(0.0, f64::MAX)) {
            assert_eq!(hit.t, 2.5);
            assert_eq!(hit.point, Vec3::new(0.0, -2.5, 0.0));
        }
    }

    #[test_log::test]
    fn test_create_traingle() {
        let m1: Material = Material::Lambertian(Lambertian::new(Color::new(1.0, 1.0, 1.0)));
        let t: Triangle = Triangle::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 2.0, 0.0), Vec3::new(2.0, 2.0, 0.0), m1);
        assert_approx_eq!(t.v0, Vec3::new(0.0, 0.0, 0.0));
        assert_approx_eq!(t.v1, Vec3::new(0.0, 2.0, 0.0));
        assert_approx_eq!(t.v2, Vec3::new(2.0, 2.0, 0.0));
    }


    #[test_log::test]
    fn test_hit_traingle() {
        let m1: Material = Material::Lambertian(Lambertian::new(Color::new(1.0, 1.0, 1.0)));
        let t: Triangle = Triangle::new(Vec3::new(-2.0, -2.0, -5.0), 
                                        Vec3::new(-2.0, 2.0, -5.0), 
                                        Vec3::new(2.0, 2.0, -5.0), 
                                        m1
                                    );

        let r: Ray = Ray::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, -1.0));
        
        match t.hit(&r, &mut Interval::new(0.0, f64::MAX)) {
            Some(hit) => {
                assert_eq!(hit.t, 5.0);
            },
            _ => { 
                panic!("Triangle should be hit")
            }
        }
    }

    #[test_log::test]
    fn test_create_sphere() {
        let m1: Material = Material::Lambertian(Lambertian::new(Color::new(1.0, 1.0, 1.0)));
        let s: Sphere = Sphere::new(Vec3::new(1.0, 2.0, -1.0), 2.0, m1);
        assert_approx_eq!(s.center, Vec3::new(1.0, 2.0, -1.0));
        assert_approx_eq!(s.radius, 2.0);
    }

    #[test_log::test]
    fn test_sphere_aabb() {
        let m1: Material = Material::Lambertian(Lambertian::new(Color::new(1.0, 1.0, 1.0)));
        let s: Sphere = Sphere::new(Vec3::new(0.0, 0.0, 0.0), 1.0, m1);

        assert_approx_eq!(s.bounding_box().axis(0).interval_min,-1.0);
        assert_approx_eq!(s.bounding_box().axis(0).interval_max,1.0);

        assert_approx_eq!(s.bounding_box().axis(1).interval_min,-1.0);
        assert_approx_eq!(s.bounding_box().axis(1).interval_max,1.0);

        assert_approx_eq!(s.bounding_box().axis(2).interval_min,-1.0);
        assert_approx_eq!(s.bounding_box().axis(2).interval_max,1.0);

        let s2: Sphere = Sphere::new(Vec3::new(2.0, 1.0, -1.0), 1.0, m1);

        assert_approx_eq!(s2.bounding_box().axis(0).interval_min,1.0);
        assert_approx_eq!(s2.bounding_box().axis(0).interval_max,3.0);

        assert_approx_eq!(s2.bounding_box().axis(1).interval_min,0.0);
        assert_approx_eq!(s2.bounding_box().axis(1).interval_max,2.0);
         
        assert_approx_eq!(s2.bounding_box().axis(2).interval_min,-2.0);
        assert_approx_eq!(s2.bounding_box().axis(2).interval_max,0.0);                   
    }

    #[test_log::test]
    fn test_hit_sphere() {
        let m1: Material = Material::Lambertian(Lambertian::new(Color::new(1.0, 1.0, 1.0)));
        let r: Ray = Ray::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, -1.0));
        let s: Sphere = Sphere::new(Vec3::new(0.0, 0.0, -3.0), 1.0, m1);

        match s.hit(&r, &mut Interval::new(0.0, f64::MAX)) {
            Some(hit) => {
                assert_eq!(hit.t, 2.0);
                assert_eq!(hit.point, Vec3::new(0.0, 0.0, -2.0));
            },
            _ => { 
                panic!("Sphere should be hit")
            }
        }
    }


}
