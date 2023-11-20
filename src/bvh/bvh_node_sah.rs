use std::f64::EPSILON;

use super::aabb::Aabb;
use crate::{
    elements::*,
    render::{Interval, Ray, Axis},
};

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
    },
}

#[allow(non_camel_case_types)]
pub struct BVH_SAH<'a> {
    /// The list of nodes of the BVH.
    pub nodes: Vec<BVHNode_SAH>,
    pub objects: &'a Vec<Element>,
}

impl BVH_SAH<'_> {
    pub fn build(objects: &Vec<Element>) -> BVH_SAH {
        // get a vec for all the indices
        let indices = (0..objects.len()).collect::<Vec<usize>>();

        // prepare the node vector
        let mut nodes = Vec::with_capacity(objects.len() * 2);

        // start with the first node, and build from there
        BVHNode_SAH::build(objects, &indices, &mut nodes, 0, 0);

        // return the tree, with a reference to the objects (needed for the hits)
        BVH_SAH {
            nodes,
            objects: &objects,
        }
    }

    pub fn hit(&self, ray: &Ray, ray_t: &mut Interval, follow: bool) -> Option<HitRecord> {
        // start at the top node, and do the hit there. this will recursively go down
        self.nodes[0].hit(self.objects, &self.nodes, ray, ray_t, follow)
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

    pub fn hit(
        &self,
        objects: &Vec<Element>,
        nodes: &Vec<BVHNode_SAH>,
        ray: &Ray,
        ray_t: &mut Interval,
        follow: bool,
    ) -> Option<HitRecord> {
        match self {
            &BVHNode_SAH::Leaf {
                parent_index: _self_parent_index,
                depth: _self_depth,
                element_index: self_element_index,
            } => {
                // we made it all the way to the leaf, so we do the hit with the actual element
                return objects[self_element_index].hit(ray, ray_t);
            }
            &BVHNode_SAH::Node {
                parent_index: _self_parent_index,
                depth: _self_depth,
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

                if self_child_l_aabb.hit(ray, ray_t) {
                    //check a hit in the following node

                    let left_hit: Option<HitRecord> =
                        nodes[self_child_l_index].hit(objects, nodes, ray, ray_t, follow);

                    match left_hit {
                        Some(lh) => {
                            ray_t.interval_max = lh.t;

                            let right_hit: Option<HitRecord> =
                                nodes[self_child_r_index].hit(objects, nodes, ray, ray_t, follow);

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
                                return nodes[self_child_r_index]
                                    .hit(objects, nodes, ray, ray_t, follow);
                            } else {
                                None
                            }
                        }
                    }
                } else {
                    // no hit on the Aabb on the left side, so we try the right
                    // first the AABB, and then a hit with the node

                    if self_child_r_aabb.hit(ray, ray_t) {
                        return nodes[self_child_r_index].hit(objects, nodes, ray, ray_t, follow);
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
            0 => 0,
            1 => {
                // only one more indice left, so we pick the first one
                let element_index = indices[0];

                // this will be the index of the new leaf
                let node_index = nodes.len();

                // add the leaf
                nodes.push(BVHNode_SAH::Leaf {
                    parent_index,
                    depth,
                    element_index,
                });

                // return the number of the node index
                return node_index;
            }
            2 => {
                // we have two elements, so get the indices of the two elements
                let element_l_indice = indices[0];
                let element_r_indice = indices[1];

                // get the index of the root node
                let node_index = nodes.len();

                // get the aabb of the roots
                let child_l_aabb = objects[element_l_indice].bounding_box();
                let child_r_aabb = objects[element_r_indice].bounding_box();

                // add the root node, and as we push the leafs next, we now the indices of those toos
                nodes.push(BVHNode_SAH::Node {
                    parent_index,
                    depth,
                    child_l_aabb,
                    child_l_index: node_index + 1,
                    child_r_aabb,
                    child_r_index: node_index + 2,
                });

                // add the first leaf
                nodes.push(BVHNode_SAH::Leaf {
                    parent_index: node_index,
                    depth: depth + 1,
                    element_index: element_l_indice,
                });

                // add the second leaf
                nodes.push(BVHNode_SAH::Leaf {
                    parent_index: node_index,
                    depth: depth + 1,
                    element_index: element_r_indice,
                });

                // return the index of the root node
                return node_index;
            }
            _ => {
                // From here on we handle the recursive case. This dummy is required, because the children
                // must know their parent, and it's easier to update one parent node than the child nodes.
                let node_index = nodes.len();
                nodes.push(BVHNode_SAH::create_dummy());

                //create the bounding box for the centroids
                let mut aabb_centroids = Aabb::new_from_points(
                    objects[indices[0]].bounding_box().centroid,
                    objects[indices[1]].bounding_box().centroid,
                );

                //loop trough all the indices, and expand the bounding box wit the the new cntroid

                for index in indices {
                    aabb_centroids = Aabb::new_from_aabbs(
                        aabb_centroids,
                        Aabb::new_from_point(objects[*index].bounding_box().centroid),
                    );
                }

                // prepare the vectors for the left and right side indicise
                let mut child_l_indices: Vec<usize> = Vec::with_capacity(indices.len() / 2);
                let mut child_r_indices: Vec<usize> = Vec::with_capacity(indices.len() / 2);

                // find the longest axis, to decide how to sort
                let sort_axis: Axis = aabb_centroids.largest_axis();
                let sort_axis_size: f64 = aabb_centroids.max[sort_axis] - aabb_centroids.min[sort_axis];

                if sort_axis_size < EPSILON {
                    // axis is too small, lets just divide into 2

                    let mid: usize = indices.len() / 2 as usize;

                    // TODO: write this one more idiomatic...
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
                    log::trace!("size: {}, splitpoint {}", sort_axis_size, splitpoint);

                    for index in indices {
                        log::trace!(
                            "point {}",
                            objects[*index].bounding_box().centroid[sort_axis]
                        );

                        if objects[*index].bounding_box().centroid[sort_axis] < splitpoint {
                            child_l_indices.push(*index);
                        } else {
                            child_r_indices.push(*index);
                        }
                    }
                }

                log::trace!(
                    "l indices {:?}, r indices {:?}",
                    child_l_indices,
                    child_r_indices
                );

                // build the left and right side nodes, and get the index back
                let child_l_index =
                    BVHNode_SAH::build(objects, &child_l_indices, nodes, node_index, depth + 1);
                let child_r_index =
                    BVHNode_SAH::build(objects, &child_r_indices, nodes, node_index, depth + 1);

                // get the bounding boxes of all the underlying objects (this we could likely do smarter)
                let child_l_aabb = Aabb::new_from_objects(objects, &child_l_indices);
                let child_r_aabb = Aabb::new_from_objects(objects, &child_r_indices);

                // push the new node
                nodes[node_index] = BVHNode_SAH::Node {
                    parent_index,
                    depth,
                    child_l_aabb,
                    child_l_index,
                    child_r_aabb,
                    child_r_index,
                };

                // return the node index
                node_index
            }
        }
    }
}
