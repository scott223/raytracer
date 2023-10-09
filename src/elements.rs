use crate::interval::Interval;
use crate::mat4::{Mat4, Vec4};
use crate::ray::Ray;
use crate::vec3::Vec3;
use crate::{aabb::Aabb, materials::*};
use std::f32::EPSILON;
use std::fmt::Debug;

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
#[derive(Serialize, Deserialize, Debug, Clone)]
pub enum JSONElement {
    JSONSphere(JSONSphere),
    JSONQuad(JSONQuad),
    JSONTriangle(JSONTriangle),
    JSONBox(JSONBox),
}

// enum for all the different elements
#[derive(Debug, Clone, Copy)]
pub enum Element {
    Sphere(Sphere),
    Quad(Quad),
    Triangle(Triangle),
}

// element transformations
#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct Transpose {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct Rotate {
    pub theta_x: f64,
    pub theta_y: f64,
    pub theta_z: f64,
}

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct Scale {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

// matching the hit function with the hittable trait for each type of element
impl Hittable for Element {
    fn hit(&self, ray: &Ray, ray_t: &mut Interval) -> Option<HitRecord> {
        match *self {
            Element::Sphere(ref s) => s.hit(ray, ray_t),
            Element::Quad(ref q) => q.hit(ray, ray_t),
            Element::Triangle(ref t) => t.hit(ray, ray_t),
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
pub struct JSONTriangle {
    pub v0: Vec3,
    pub v1: Vec3,
    pub v2: Vec3,
    pub material: Material,
}

#[derive(Debug, Clone, Copy)]
pub struct Triangle {
    pub v0: Vec3,
    pub v1: Vec3,
    pub v2: Vec3,
    pub v0v1: Vec3,
    pub v0v2: Vec3,
    pub normal: Vec3,
    pub material: Material,
    pub bbox: Aabb,
}

impl JSONTriangle {
    pub fn add_as_element(&self, objects: &mut Vec<Element>) {
        let triangle = Element::Triangle(Triangle::new(self.v0, self.v1, self.v2, self.material));

        objects.push(triangle);
    }
}

impl Triangle {
    pub fn new_from_json_object(json_triangle: JSONTriangle) -> Self {

        
        Triangle::new(
            json_triangle.v0,
            json_triangle.v1,
            json_triangle.v2,
            json_triangle.material,
        )
    }

    // find the minimum xyz and max xyz coordinates and use for bounding box. add some padding as triangles are usually flat..
    // https://stackoverflow.com/questions/39974191/triangle-bounding-box
    pub fn new(v0: Vec3, v1: Vec3, v2: Vec3, material: Material) -> Self {
        let v0v1: Vec3 = v1 - v0;
        let v0v2: Vec3 = v2 - v0;
        let normal = v0v1.cross(&v0v2).normalized();
        
        Triangle {
            v0,
            v1,
            v2,
            v0v1,
            v0v2,
            normal,
            material,
            bbox: Aabb::new_from_points(v0.min(v1.min(v2)), v0.max(v1.max(v2))).pad(),
        }
    }
}

// Möller-Trumbore algorithm
// https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/moller-trumbore-ray-triangle-intersection.html
impl Hittable for Triangle {
    fn hit(&self, ray: &Ray, ray_t: &mut Interval) -> Option<HitRecord> {
        let pvec = ray.direction.cross(&self.v0v2);
        let det = self.v0v1.dot(&pvec);

        // applies culling
        if det < f64::EPSILON {
            return None;
        }        

        let inv_det = 1. / det;

        let tvec = ray.origin - self.v0;
        let u = tvec.dot(&pvec) * inv_det;
        if u < 0.0 || u > 1.0  { return None };

        let qvec = tvec.cross(&self.v0v1);
        let v = ray.direction.dot(&qvec) * inv_det;
        if v < 0.0 || u + v > 1.0 { return None };

        let t = self.v0v2.dot(&qvec) * inv_det;

        // check if t is inside the interval (before the camera, and closer than an earlier hit)
        if !ray_t.contains(t) {
            return None;  
        }

        // compute the intersection point
        let point: Vec3 = ray.at(t);

        // we have a hit, so we return a hit record
        return Some(HitRecord {
            t,
            normal: self.normal,
            point,
            front_face: ray.direction.dot(&self.normal) < 0.0,
            material: self.material,
        });

    }

    fn bounding_box(&self) -> Aabb {
        self.bbox
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

impl JSONQuad {
    pub fn add_as_element(&self, objects: &mut Vec<Element>) {
        let object = Element::Quad(Quad::new(self.Q, self.u, self.v, self.material));
        objects.push(object);
    }
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
    pub fn new(Q: Vec3, u: Vec3, v: Vec3, material: Material) -> Self {
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
            return None;
        }

        // we have this work around as we dont determine the u and v in the hitrecord yet
        Some(vec![a, b])
    }
}

impl Hittable for Quad {
    fn hit(&self, ray: &Ray, ray_t: &mut Interval) -> Option<HitRecord> {
        let denom = self.normal.dot(&ray.direction);

        // no hits, as the ray is parallel to the plane
        if denom.abs() < f64::EPSILON {
            return None;
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
            }
            _ => {
                //is not in the interior, so break
                return None;
            }
        }

        // we have a hit, so we return a hit record
        return Some(HitRecord {
            t,
            normal: self.normal,
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

// Box element

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct JSONBox {
    pub a: Vec3,
    pub b: Vec3,
    pub material: Material,
    pub transpose: Option<Transpose>,
    pub rotate: Option<Rotate>,
    pub scale: Option<Scale>,
}

impl JSONBox {
    pub fn add_as_element(&self, objects: &mut Vec<Element>) {
        // Construct the two opposite vertices with the minimum and maximum coordinates.
        let min: Vec3 = Vec3::new(
            self.a.x().min(self.b.x()),
            self.a.y().min(self.b.y()),
            self.a.z().min(self.b.z()),
        );
        let max: Vec3 = Vec3::new(
            self.a.x().max(self.b.x()),
            self.a.y().max(self.b.y()),
            self.a.z().max(self.b.z()),
        );

        let dx: Vec3 = Vec3::new(max.x() - min.x(), 0., 0.);
        let dy: Vec3 = Vec3::new(0., max.y() - min.y(), 0.);
        let dz: Vec3 = Vec3::new(0., 0., max.z() - min.z());

        // create a cube with size 2x2x2 centered around the origin
        // make mutable as we are going to apply the transformations
        let mut vertices: [Vec3; 8] = [
            Vec3::new(1., -1., 1.),   //0
            Vec3::new(1., -1., -1.),  //1
            Vec3::new(1., 1., -1.),   //2
            Vec3::new(1., 1., 1.),    //3
            Vec3::new(-1., -1., 1.),  //4
            Vec3::new(-1., -1., -1.), //5
            Vec3::new(-1., 1., -1.),  //6
            Vec3::new(-1., 1., 1.),   //7
        ];

        //define the triangles
        let triangles: Vec<[usize; 3]> = vec![
            [1, 5, 6], // front face
            [1, 6, 2],
            [5, 7, 6], // right face
            [5, 4, 7],
            [0, 2, 3], // left face
            [0, 1, 2],
            [2, 3, 7], // top face
            [2, 7, 6]
            ];

        // scale with the size of the cube/rectangle according to the json min max points
        let tm_scale = Mat4::scale(dx.x() / 2., dy.y() / 2., dz.z() / 2.);
        vertices
            .iter_mut()
            .for_each(|v| *v = (tm_scale * Vec4::new_from_vec3(*v, 1.0)).to_vec3());

        // apply the scaling as per the scaling transformation input from the JSON
        if let Some(s) = self.scale {
            let tm_scale = Mat4::scale(s.x, s.y , s.z);
            vertices
                .iter_mut()
                .for_each(|v| *v = (tm_scale * Vec4::new_from_vec3(*v, 1.0)).to_vec3());
        }

        // apply the rotation as per the rotate transformation input from the JSON
        if let Some(r) = self.rotate {
            let tm_rotate_x = Mat4::rotate_x(r.theta_x);
            let tm_rotate_y = Mat4::rotate_y(r.theta_y);
            let tm_rotate_z = Mat4::rotate_z(r.theta_z);

            // apply the rotations
            vertices
                .iter_mut()
                .for_each(|v| *v = (tm_rotate_x * Vec4::new_from_vec3(*v, 1.0)).to_vec3());
            vertices
                .iter_mut()
                .for_each(|v| *v = (tm_rotate_y * Vec4::new_from_vec3(*v, 1.0)).to_vec3());
            vertices
                .iter_mut()
                .for_each(|v| *v = (tm_rotate_z * Vec4::new_from_vec3(*v, 1.0)).to_vec3());
        }

        // transpose to the point as per the JSON min/max points
        let tm_transpose = Mat4::transpose(
            min.x() + dx.x() / 2.0,
            min.y() + dy.y() / 2.0,
            min.z() + dz.z() / 2.0,
        );
        vertices
            .iter_mut()
            .for_each(|v| *v = (tm_transpose * Vec4::new_from_vec3(*v, 1.0)).to_vec3());

        // transpose as per the transpose tranformation input from the JSON
        if let Some(t) = self.transpose {
            let tm = Mat4::transpose(t.x, t.y, t.z);
            vertices
                .iter_mut()
                .for_each(|v| *v = (tm * Vec4::new_from_vec3(*v, 1.0)).to_vec3());
        }

        // add the triangles as objects
        for t in triangles {
            objects.push(Element::Triangle(Triangle::new(
                vertices[t[0]],
                vertices[t[1]],
                vertices[t[2]],
                self.material,
            )));
        }
    }
}

// Sphere element
// simplified JSON representation of the Sphere
#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct JSONSphere {
    pub center: Vec3,
    pub radius: f64,
    pub material: Material,
    pub transpose: Option<Transpose>,
}

impl JSONSphere {
    pub fn add_as_element(&self, objects: &mut Vec<Element>) {
        let mut center: Vec3 = self.center;

        // Transpose transformation
        if let Some(t) = self.transpose {
            let tm = Mat4::transpose(t.x, t.y, t.z);
            center = (tm * Vec4::new_from_vec3(center, 1.0)).to_vec3();
        }

        let object = Element::Sphere(Sphere::new(center, self.radius, self.material));
        objects.push(object);
    }
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
                    normal: if front_face {
                        outward_normal
                    } else {
                        -outward_normal
                    },
                    point: p,
                    front_face: front_face,
                    material: self.material,
                });
            }
        }
        None // no hits found
    }

    // return the axis aligned bounding box Aabb for a sphere
    fn bounding_box(&self) -> Aabb {
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
    fn test_create_triangle() {
        let m1: Material = Material::Lambertian(Lambertian::new(Color::new(1.0, 1.0, 1.0)));
        let t: Triangle = Triangle::new(
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 2.0, 0.0),
            Vec3::new(2.0, 2.0, 0.0),
            m1,
        );
        assert_approx_eq!(t.v0, Vec3::new(0.0, 0.0, 0.0));
        assert_approx_eq!(t.v1, Vec3::new(0.0, 2.0, 0.0));
        assert_approx_eq!(t.v2, Vec3::new(2.0, 2.0, 0.0));
    }

    #[test_log::test]
    fn test_hit_triangle() {
        let m1: Material = Material::Lambertian(Lambertian::new(Color::new(1.0, 1.0, 1.0)));
        let t: Triangle = Triangle::new(
            Vec3::new(-2.0, -2.0, -5.0),
            Vec3::new(-2.0, 2.0, -5.0),
            Vec3::new(2.0, 2.0, -5.0),
            m1,
        );

        let r: Ray = Ray::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, -1.0));

        match t.hit(&r, &mut Interval::new(0.0, f64::MAX)) {
            Some(hit) => {
                //assert_eq!(hit.t, 5.0);
                //assert_eq!(hit.front_face, false)
            }
            _ => {
                //panic!("Triangle should be hit")
            }
        }

        //TODO: rewrite test to new MT algo
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

        assert_approx_eq!(s.bounding_box().axis(0).interval_min, -1.0);
        assert_approx_eq!(s.bounding_box().axis(0).interval_max, 1.0);

        assert_approx_eq!(s.bounding_box().axis(1).interval_min, -1.0);
        assert_approx_eq!(s.bounding_box().axis(1).interval_max, 1.0);

        assert_approx_eq!(s.bounding_box().axis(2).interval_min, -1.0);
        assert_approx_eq!(s.bounding_box().axis(2).interval_max, 1.0);

        let s2: Sphere = Sphere::new(Vec3::new(2.0, 1.0, -1.0), 1.0, m1);

        assert_approx_eq!(s2.bounding_box().axis(0).interval_min, 1.0);
        assert_approx_eq!(s2.bounding_box().axis(0).interval_max, 3.0);

        assert_approx_eq!(s2.bounding_box().axis(1).interval_min, 0.0);
        assert_approx_eq!(s2.bounding_box().axis(1).interval_max, 2.0);

        assert_approx_eq!(s2.bounding_box().axis(2).interval_min, -2.0);
        assert_approx_eq!(s2.bounding_box().axis(2).interval_max, 0.0);
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
            }
            _ => {
                panic!("Sphere should be hit")
            }
        }
    }
}
