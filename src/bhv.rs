use std::fmt::Debug;

use crate::aabb::Aabb;
use crate::elements::*;
use crate::interval::Interval;
use crate::ray::Ray;
use crate::vec3::Vec3;

use rand::{Rng, SeedableRng};
use rand::rngs::SmallRng;

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
        //log::info!("start: {}, end: {}, span: {}", start, end, object_span);

        if object_span == 1 {
            // we just have one element, so we can add that as a final node
            // log::info!(
            //    "one element left, returning that element as a final node. object n: {}",
            //    start
            //);
            return Box::new(objects[start]);
        } else if object_span == 2 {
            // we have two items, lets see which one comes first and asign in the right order to the end node
            if BHVNode::box_compare(&objects[start], &objects[start + 1]) {
                let node: BHVNode = BHVNode {
                    left: Box::new(objects[start]),
                    right: Box::new(objects[start + 1]),
                    bbox: Aabb::new_from_aabbs(
                        objects[start].bounding_box(),
                        objects[start + 1].bounding_box(),
                    ),
                };

                // log::info!(
                //    "Two items, creating end node (BHV with two elements) for n: {} and n+1: {}",
                //    start,
                //    start + 1
                //);
                return Box::new(node);
            } else {
                let node: BHVNode = BHVNode {
                    left: Box::new(objects[start + 1]),
                    right: Box::new(objects[start]),
                    bbox: Aabb::new_from_aabbs(
                        objects[start + 1].bounding_box(),
                        objects[start].bounding_box(),
                    ),
                };

                // log::info!("Two items, flipping & creating end node (BHV with two elements) for n: {} and n+1: {}", start+1, start);
                return Box::new(node);
            }
        } else {
            // we still have a few elements, so we need to create some more nodes and pass the elements
            // we randomize the axis that we sort on each time
            let n = SmallRng::seed_from_u64(223).gen_range(0..3);
            // log::info!("Sorting over {} axis", n);
            objects[start..end].sort_by(|a, b| {
                a.bounding_box()
                    .axis(n)
                    .interval_min
                    .total_cmp(&b.bounding_box().axis(n).interval_min)
            });

            // we cut the sample in half
            let mid: usize = start + object_span / 2;

            // we create two sub nodes, that get the same object list, but a will look a at a smaller subsection
            let left = BHVNode::new(objects, start, mid);
            let right = BHVNode::new(objects, mid, end);

            // get the bounding box
            let bbox = Aabb::new_from_aabbs(left.bounding_box(), right.bounding_box());

            // create the node
            let node: BHVNode = BHVNode {
                left,
                right,
                bbox,
            };

            // log::info!(
            //    "More than two items (s: {}, m: {}, e: {}), creating a new layer of BHV nodes.",
            //    start,
            //    mid,
            //    end
            //);
            return Box::new(node);
        }
    }

    // compare function to sort, currently only breaks down along the y-axis (this is axis n=1)
    pub fn box_compare(a: &dyn Hittable, b: &dyn Hittable) -> bool {
        return a.bounding_box().axis(1).interval_min < b.bounding_box().axis(1).interval_min;
    }
}

impl Debug for BHVNode {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "left: {:?}, right: {:?}, aabb: {:?}",
            self.left,
            self.right,
            self.bounding_box()
        )
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
        // this is where the real speadup happens, as we dont have to do any fancy calculations, just check for big box at the start of the tree (of further down below)
        if !self.bounding_box().hit(ray, ray_t) {
            return None;
        }

        // check for a hit in the left path first
        let left_hit = self.left.hit(ray, ray_t);

        //TODO we can do the next part smoother, where we immediatly do the right hit too (adjust the interval before) and then just make one match statement

        match left_hit {
            Some(lh) => {
                //there is a hit on the left path, so lets adjust the interval and see if we have one closer by on the right
                ray_t.interval_max = lh.t;
                let right_hit = self.right.hit(ray, ray_t);

                match right_hit {
                    Some(rh) => {
                        // we have a closer hit on the right, so return the right hit
                        return Some(rh);
                    }
                    _ => {
                        // there is no closer hit on the right, so return the left hit
                        return Some(lh);
                    }
                }
            }
            _ => {
                // no hit on the left side, so lets try the right with the unmodified interval ray_t
                let right_hit = self.right.hit(ray, ray_t);
                match right_hit {
                    Some(rh) => {
                        // there is a hit on the right, so we return that hit
                        return Some(rh);
                    }
                    _ => {
                        // no hit on the left or right, so we return a Nones
                        return None;
                    }
                }
            }
        }
    }

    fn bounding_box(&self) -> Aabb {
        self.bbox //return the pre-set bounding box
    }

    // need to implement these for the Hittable trait, but serve little function in the BHVNode
    fn pdf_value(&self, _origin: Vec3, _direction: Vec3) -> f64 {
        0.0
    }
    
    // need to implement these for the Hittable trait, but serve little function in the BHVNode
    fn random(&self, _origin: Vec3) -> Vec3 {
        Vec3::new(1.0, 0.0, 0.0)
    }
}

//TODO: Tests