use std::fmt::Debug;
use crate::{materials::*, aabb::Aabb};
use crate::ray::Ray;
use crate::vec3::Vec3;
use crate::interval::Interval;

use rand::Rng;
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

// enum for all the different elements
#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub enum Element {
    Sphere(Sphere),
//    Plane(Plane),
}

// matching the hit function with the hittable trait for each type of element
impl Hittable for Element {
    fn hit(&self, ray: &Ray, ray_t: &mut Interval) -> Option<HitRecord> {
        match *self {
            Element::Sphere(ref s) => s.hit(ray, ray_t),
 //           Element::Plane(ref p) => p.hit(ray, t_min, t_max),
        }
    }

    fn bounding_box(&self) -> Aabb {
        match *self {
            Element::Sphere(ref s) => s.bounding_box(),
        }
    }
}

// this is a node for the BHV tree, and it consists of objects with a hittable trait (either another node, or an element)
// it will only have a left or a right object
pub struct BHVNode {
    left: Box<dyn Hittable + Sync>,
    right: Box<dyn Hittable + Sync>,
    bbox: Aabb, // this gets calculated and set when constructing the tree, so we can use this to speed up the hit calculations
}

impl BHVNode {
    pub fn new(objects: &mut Vec<Element>, start: usize, end: usize) -> Box<dyn Hittable + Sync> {

        // see how many elements we still have left
        let object_span: usize = end - start;
        log::info!("start: {}, end: {}, span: {}", start, end, object_span);
        
        if object_span == 1 {
            // we just have one element, so we can add that as a final node
            println!("one element left, returning that element as a final node. object n: {}", start);
            return Box::new(objects[start]);

        } else if object_span == 2 {
            // we have two items, lets see which one comes first and asign in the right order to the end node
            if BHVNode::box_compare(&objects[start], &objects[start+1]) {

                let node: BHVNode = BHVNode {
                    left: Box::new(objects[start]),
                    right: Box::new(objects[start+1]),
                    bbox: Aabb::new_from_aabbs(objects[start].bounding_box(), objects[start+1].bounding_box()),
                };

                log::info!("Two items, creating end node (BHV with two elements) for n: {} and n+1: {}", start, start+1);
                return Box::new(node);
            } else {
                let node: BHVNode = BHVNode {
                    left: Box::new(objects[start+1]),
                    right: Box::new(objects[start]),
                    bbox: Aabb::new_from_aabbs(objects[start+1].bounding_box(), objects[start].bounding_box()),
                };

                log::info!("Two items, flipping & creating end node (BHV with two elements) for n: {} and n+1: {}", start+1, start);
                return Box::new(node);               
            }
            
        } else {
            // we still have a few elements, so we need to create some more nodes and pass the elements
            // we randomize the axis that we sort on each time
            let n = rand::thread_rng().gen_range(0..3);
            log::info!("Sorting over {} axis", n);
            objects[start..end].sort_by(|a, b| a.bounding_box().axis(n).interval_min.total_cmp(&b.bounding_box().axis(n).interval_min));

            // we cut the sample in half
            let mid: usize = start + object_span/2;

            // we create two sub nodes, that get the same object list, but a will look a at a smaller subsection
            let left = BHVNode::new(objects, start, mid);
            let right = BHVNode::new(objects, mid, end);

            // get the bounding box
            let bbox = Aabb::new_from_aabbs(left.bounding_box(), right.bounding_box());

            // create the node
            let node: BHVNode = BHVNode {
                left: left,
                right: right,
                bbox,
            };

            log::info!("More than two items (s: {}, m: {}, e: {}), creating a new layer of BHV nodes.", start, mid, end);
            return Box::new(node);

        }
    }

    // compare function to sort, currently only breaks down along the z-axis (this is axis n=2)
    pub fn box_compare(a: &dyn Hittable, b: &dyn Hittable) -> bool {
        return a.bounding_box().axis(2).interval_min < b.bounding_box().axis(2).interval_min
    }

}

// display trait
//impl fmt::Display for BHVNode {
//    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
//      write!(f, "[left: {:?}, right: {:?}, aabb: {:?}]", self.left, self.right, self.aabb)
//    }
//}

impl Debug for BHVNode {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "left: {:?}, right: {:?}, aabb: {:?}", self.left, self.right, self.bounding_box())
    }
}

impl Debug for dyn Hittable + Sync {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "aabb: {:?}", self.bounding_box())
    }
}

impl Hittable for BHVNode {

    // recursive check for hits through the BHV nodes
    fn hit(&self, ray: &Ray, ray_t: &mut Interval) -> Option<HitRecord> {
        // check if we hit the bounding box of this node, because if we dont, we can stop right here
        // this is where the real speadup happens, as we dont have to do any fancy calculations, just check for big box
        if !self.bounding_box().hit(ray, ray_t) { 
           return None 
        }

        // check for a hit in the left path first
        let left_hit = self.left.hit(ray, ray_t);
        
        match left_hit {
            Some(lh) => {
                //there is a hit on the left path, so lets adjust the interval and see if we have one closer by on the right
                ray_t.interval_max = lh.t;
                let right_hit = self.right.hit(ray, ray_t);

                match right_hit {
                    Some(rh) => {
                        // we have a closer hit on the right, so return the right hit
                        return Some(rh);
                    },
                    _ => {
                        // there is no closer hit on the right, so return the left hit
                        return Some(lh);
                    }
                }
            },
            _ => {
                // no hit on the left side, so lets try the right with the unmodified interval ray_t
                let right_hit = self.right.hit(ray, ray_t);
                match right_hit {
                    Some(rh) => {
                        // there is a hit on the right, so we return that one
                        return Some(rh);
                    },
                    _ => {
                        // no hit on the left or right, so we return nothing
                        return None;
                    }
                }
            }
        }

        //not sure if this one ever gets called, but just in case we return nothing if there is nothing
        //None

    }
    
    fn bounding_box(&self) -> Aabb {
        self.bbox
        //Aabb::new_from_aabbs(self.left.bounding_box(), self.right.bounding_box())
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

            if distance >= 0.0 {
                if distance > ray_t.interval_min && distance < ray_t.interval_max {
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
        }
        None // no hits found
    }

    fn bounding_box(&self) -> Aabb {
        // TODO, will need to figure out how to add a bounding box for an infinite plane
        todo!();
    }
}

// Sphere element

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct Sphere {
    pub center: Vec3,
    pub radius: f64,
    pub material: Material,
}

// Creating a new Sphere with a center and a radius
impl Sphere {
    pub fn new(center: Vec3, radius: f64, material: Material) -> Self {
        Sphere {
            center: center,
            radius: radius,
            material: material,
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
        let rvec = Vec3::new(self.radius, self.radius, self.radius);
        Aabb::new_from_points(self.center - rvec, self.center + rvec)
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

        if let Some(hit) = s.hit(&r, &mut Interval::new(0.0, f64::MAX)) {
            assert_eq!(hit.t, 2.0);
            assert_eq!(hit.point, Vec3::new(0.0, 0.0, -2.0));
        }
    }
}
