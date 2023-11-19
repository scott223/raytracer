use serde::Deserialize;
use serde::Serialize;

use crate::linalg::{Mat4, Vec3, Vec4};
use crate::materials::Material;

use super::Element;
use super::Rotate;
use super::Scale;
use super::Transpose;
use super::Triangle;

// Box element

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct JSONCube {
    pub a: Vec3,
    pub b: Vec3,
    pub material: Material,
    pub transpose: Option<Transpose>,
    pub rotate: Option<Rotate>,
    pub scale: Option<Scale>,
}

impl JSONCube {
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

        //         3     7
        //    2    x - x
        //      x - x6|
        //      |0x-| x 4
        //      x - x
        //    1      5
        //

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
            [2, 7, 3], // top face
            [2, 6, 7],
            [0, 7, 4], // back face
            [0, 3, 7],
            [4, 5, 1], // bottom face
            [1, 0, 4],
        ];

        // TODO: create one big transformation matrix and apply once

        // scale with the size of the cube/rectangle according to the json min max points
        let tm_scale = Mat4::scale(dx.x() / 2., dy.y() / 2., dz.z() / 2.);
        vertices
            .iter_mut()
            .for_each(|v| *v = (tm_scale * Vec4::new_from_vec3(*v, 1.0)).to_vec3());

        // apply the scaling as per the scaling transformation input from the JSON
        if let Some(s) = self.scale {
            let tm_scale = Mat4::scale(s.x, s.y, s.z);
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
            //print!("v0, v1, v2: {:.2}, {:.2}, {:.2}", vertices[t[0]], vertices[t[1]], vertices[t[2]]);
            objects.push(Element::Triangle(Triangle::new(
                vertices[t[0]],
                vertices[t[1]],
                vertices[t[2]],
                self.material,
            )));
        }
    }
}
