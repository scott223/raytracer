use std::fs::File;
use std::io::Read;

use serde::Deserialize;
use serde::Serialize;
use wavefront_obj::obj;

use crate::elements::Triangle;
use crate::linalg::{Mat4, Vec3, Vec4};
use crate::materials::Material;

use super::Element;
use super::Rotate;
use super::Scale;
use super::Transpose;

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct JSONObj {
    pub filepath: String,
    pub transpose: Option<Transpose>,
    pub rotate: Option<Rotate>,
    pub scale: Option<Scale>,
    pub material: Material,
}

impl JSONObj {
    pub fn add_as_element(&self, objects: &mut Vec<Element>) {
        //todo ERROR
        //let input = BufReader::new(File::open(self.filepath).unwrap());

        let mut f = File::open(self.filepath.to_owned()).unwrap();
        let mut s = String::new();
        let _ = f.read_to_string(&mut s);

        let result = obj::parse(s).unwrap();

        log::trace!("loaded .obj with {} objects", result.objects.len());

        for object in result.objects {
            let mut vertices: Vec<Vec3> = Vec::new();
            object
                .vertices
                .iter()
                .for_each(|v| vertices.push(Vec3::new(v.x, v.y, v.z)));

            //TODO move the scaling, rotation and transpose to a seperate functio
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

            if let Some(t) = self.transpose {
                let tm = Mat4::transpose(t.x, t.y, t.z);
                vertices
                    .iter_mut()
                    .for_each(|v| *v = (tm * Vec4::new_from_vec3(*v, 1.0)).to_vec3());
            }

            for geometry in object.geometry {
                log::trace!("Loaded {} shapes", geometry.shapes.len());
                for shape in geometry.shapes {
                    match shape.primitive {
                        wavefront_obj::obj::Primitive::Triangle(a, b, c) => {
                            let triangle = Element::Triangle(Triangle::new(
                                Vec3::new(vertices[a.0].x(), vertices[a.0].y(), vertices[a.0].z()),
                                Vec3::new(vertices[b.0].x(), vertices[b.0].y(), vertices[b.0].z()),
                                Vec3::new(vertices[c.0].x(), vertices[c.0].y(), vertices[c.0].z()),
                                self.material,
                            ));

                            objects.push(triangle);
                        }
                        _ => {}
                    }
                }
            }
        }

        //objects.push(Element::Triangle(Triangle::new(
        // vertices[t[0]],
        // vertices[t[1]],
        //   vertices[t[2]],
        //     self.material,
        //   )));
    }
}
