

use std::f64::EPSILON;

use super::aabb::Aabb;
use crate::{elements::*, render::{Interval, Ray}};


// this is a node for the BHV tree, and it consists of objects with a hittable trait (either another node, or an element)
// it will only have a left or a right object
#[allow(non_camel_case_types)]

pub enum BVHNode_SAH {
    Leaf {
        parent_index: usize,
        depth: u32,
        element_index: usize,
    },
    Node {
        parent_index: usize,
        depth: u32,
        /// Index of the left subtree's root node.
        child_l_index: usize,

        /// The aabb of the left element.
        child_l_aabb: Aabb,

        /// Index of the right subtree's root node.
        child_r_index: usize,

        /// The aabb of the left element.
        child_r_aabb: Aabb,
    }
}

pub struct BVH_SAH {
    /// The list of nodes of the [`BVH`].
    ///
    /// [`BVH`]: struct.BVH.html
    ///
    pub nodes: Vec<BVHNode_SAH>,
}

impl BVH_SAH {
    pub fn build(objects: &Vec<Element>) -> BVH_SAH {
        let indices = (0..objects.len()).collect::<Vec<usize>>();
        let expected_node_count = objects.len() * 2;
        let mut nodes = Vec::with_capacity(expected_node_count);
        BVHNode_SAH::build(objects, &indices, &mut nodes, 0, 0);
        BVH_SAH { nodes }
    }

}

impl BVHNode_SAH {

    /// The build function sometimes needs to add nodes while their data is not available yet.
    /// A dummy created by this function serves the purpose of being changed later on.
    fn create_dummy() -> BVHNode_SAH {
        BVHNode_SAH::Leaf {
            parent_index: 0,
            depth: 0,
            element_index: 0,
        }
    }

    pub fn hit(&self, objects: &Vec<Element>, nodes: &Vec<BVHNode_SAH>, ray: &Ray, ray_t: &mut Interval) -> Option<HitRecord> {

        match (self) {
            &BVHNode_SAH::Leaf {
                parent_index: self_parent_index,
                depth: self_depth,
                element_index: self_element_index,
            } => {
                
                // we made it all the way to the leaf, so we do the hit with the actual object
                return objects[self_element_index].hit(ray, ray_t);

            },
            &BVHNode_SAH::Node{
                parent_index: self_parent_index,
                depth: self_depth,
                /// Index of the left subtree's root node.
                child_l_index: self_child_l_index,
        
                /// The aabb of the left element.
                child_l_aabb: self_child_l_aabb,
        
                /// Index of the right subtree's root node.
                child_r_index: self_child_r_index,
        
                /// The aabb of the left element.
                child_r_aabb: self_child_r_aabb,
            } => {

                // check for a hit in the left aabb first
                //if self_child_l_aabb.hit(ray, ray_t) {

                if self_child_l_aabb.hit(ray, ray_t) {

                    //check a hit in the following node

                    let left_hit: Option<HitRecord> = nodes[self_child_l_index].hit(objects, nodes, ray, ray_t);

                    match left_hit {
                        Some(lh) => {

                            ray_t.interval_max = lh.t;
                            let right_hit: Option<HitRecord> = nodes[self_child_r_index].hit(objects, nodes, ray, ray_t);

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
                            // this function returns Some(rh) if a hit is found, and None if no hit is found
                            
                            if self_child_r_aabb.hit(ray, ray_t) {
                                return nodes[self_child_r_index].hit(objects, nodes, ray, ray_t);
                            } else {
                                None
                            }
                        }
                    }

                } else {

                    if (self_child_r_aabb.hit(ray, ray_t)) {
                        return nodes[self_child_r_index].hit(objects, nodes, ray, ray_t);
                    } else {
                        None
                    }   
                }
            }
        }
    }

    pub fn build(
        objects: &Vec<Element>,
        indices: &[usize],
        nodes: &mut Vec<BVHNode_SAH>,
        parent_index: usize,
        depth: u32,
    ) -> usize {

        match indices.len() {
            0 => { 0 },
            1 => {
                // only one more indice left
                let element_index = indices[0];

                // this will be the index of the new leaf
                let node_index = nodes.len();

                // add the leaf
                nodes.push(BVHNode_SAH::Leaf {
                    parent_index,
                    depth,
                    element_index,
                });

                return node_index;
            },
            2 => {
                let element_l_indice = indices[0];
                let element_r_indice = indices[1];

                let node_index = nodes.len();

                // add the root
                let child_l_aabb = Aabb::default();
                let child_r_aabb = Aabb::default();

                nodes.push(BVHNode_SAH::Node {
                    parent_index,
                    depth,
                    child_l_aabb,
                    child_l_index: node_index+1,
                    child_r_aabb,
                    child_r_index: node_index+2,
                });

                // add the first leaf
                nodes.push(BVHNode_SAH::Leaf {
                    parent_index: node_index,
                    depth,
                    element_index: element_l_indice,
                });

                // add the second leaf
                nodes.push(BVHNode_SAH::Leaf {
                    parent_index: node_index,
                    depth,
                    element_index: element_r_indice,
                });

                return node_index;

            }
            _ => {    

                // From here on we handle the recursive case. This dummy is required, because the children
                // must know their parent, and it's easier to update one parent node than the child nodes.
                let node_index = nodes.len();
                nodes.push(BVHNode_SAH::create_dummy());

                // create the bounding box for this entire set
                let mut aabb_objects = objects[indices[0]].bounding_box();

                for index in indices {
                    aabb_objects = Aabb::new_from_aabbs(aabb_objects, objects[*index].bounding_box());
                }

                //create the bounding box for the centroids
                let mut aabb_centroids = Aabb::new_from_points(objects[indices[0]].bounding_box().centroid, objects[indices[1]].bounding_box().centroid);

                for index in indices {
                    
                    aabb_centroids = Aabb::new_from_aabbs(aabb_centroids, Aabb::new_from_points(objects[*index].bounding_box().centroid,objects[*index].bounding_box().centroid));
                }

                let mut child_l_indices: Vec<usize> = Vec::with_capacity(indices.len()/2);
                let mut child_r_indices: Vec<usize> = Vec::with_capacity(indices.len()/2);
                
                let sort_axis = aabb_centroids.largest_axis();
                let sort_axis_size = aabb_centroids.max[sort_axis] - aabb_centroids.min[sort_axis];

                if sort_axis_size < EPSILON {
                    // axis is too small, lets just divide into 2

                    let mid = indices.len() / 2 as usize;

                    for (i, index) in indices.iter().enumerate() {
                        if i < mid {
                            child_l_indices.push(*index);
                        } else {
                            child_r_indices.push(*index);
                        }
                    }

                } else {

                    //axis is big enough to split
                    // TODO here we need to add the SAH logic!

                    let splitpoint = aabb_centroids.min[sort_axis] + sort_axis_size / 2.0;
                    print!("size: {}, splitpoint {}", sort_axis_size, splitpoint);

                    for index in indices {
                        
                        println!("point {}", objects[*index].bounding_box().centroid[sort_axis]);
                        
                        if objects[*index].bounding_box().centroid[sort_axis] < splitpoint {
                            child_l_indices.push(*index);
                        } else {
                            child_r_indices.push(*index);
                        }
                        
                    }

                }   

                println!("l indices {:?}, r indices {:?}", child_l_indices, child_r_indices);         

                let child_l_index =
                BVHNode_SAH::build(objects, &child_l_indices, nodes, node_index, depth + 1);
                let child_r_index =
                BVHNode_SAH::build(objects, &child_r_indices, nodes, node_index, depth + 1);

                let mut child_l_aabb = Aabb::default();
                let mut child_r_aabb = Aabb::default();

                nodes[node_index] = BVHNode_SAH::Node {
                    parent_index,
                    depth,
                    child_l_aabb,
                    child_l_index,
                    child_r_aabb,
                    child_r_index,
                };
        
                node_index

            }
        }
    }
}

/* 
impl Hittable for BVHNode_SAH<'_> {
    // recursive check for hits through the BHV nodes
    fn hit(&self, ray: &Ray, ray_t: &mut Interval) -> Option<HitRecord> {
        // check if we hit the bounding box of this node, because if we dont, we can stop right here
        // this is where the real speadup happens, as we dont have to do any fancy calculations, just check for big box at the start of the tree (of further down below)
        if !self.bounding_box().hit(ray, ray_t) {
            return None;
        }

        let left_hit = self.left?.hit(ray, ray_t);

        //TODO we can do the next part smoother, where we immediatly do the right hit too (adjust the interval before) and then just make one match statement

        match left_hit {
            Some(lh) => {
                //there is a hit on the left path, so lets adjust the interval and see if we have one closer by on the right
                ray_t.interval_max = lh.t;
                let right_hit = self.right?.hit(ray, ray_t);

                match right_hit {
                    Some(rh) => {
                        // we have a closer hit on the right, so return the right hit
                        Some(rh)
                    }
                    _ => {
                        // there is no closer hit on the right, so return the left hit
                        Some(lh)
                    }
                }
            }
            _ => {
                // no hit on the left side, so lets try the right with the unmodified interval ray_t
                // this function returns Some(rh) if a hit is found, and None if no hit is found
                self.right?.hit(ray, ray_t)
            }
        }

    }

    fn bounding_box(&self) -> Aabb {
        self.bounding_box //return the pre-set bounding box
    }

    // need to implement these for the Hittable trait, but serve little function in the BHVNode
    fn pdf_value(&self, _origin: Vec3, _direction: Vec3) -> f64 {
        0.0
    }

    // need to implement these for the Hittable trait, but serve little function in the BHVNode
    fn random(&self, _origin: Vec3, _rng: &mut SmallRng) -> Vec3 {
        Vec3::new(1.0, 0.0, 0.0)
    }

}
*/
