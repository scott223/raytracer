use criterion::{criterion_group, criterion_main, Criterion};
use rand::seq::SliceRandom;
use rand_distr::Distribution;
use rand_distr::Uniform;
use raytracer::bvh::Aabb;

use raytracer::bvh::BVHSplitMethod;
use raytracer::elements::JSONElement;
use raytracer::elements::JSONObj;
use raytracer::elements::Rotate;
use raytracer::elements::Scale;
use raytracer::elements::Transpose;
use raytracer::linalg::Vec3;
use raytracer::materials::Lambertian;
use raytracer::render::Color;
use raytracer::render::Config;
use raytracer::render::Interval;
use raytracer::render::JSONCamera;
use raytracer::render::JSONScene;
use raytracer::render::Ray;
use raytracer::render::RenderIntegrator;

pub fn criterion_benchmark(c: &mut Criterion) {
    let bbox: Aabb = Aabb::new_from_points(Vec3::new(1.0, 1.0, 1.0), Vec3::new(2.0, 2.0, 2.0));

    let step = Uniform::new(0.0, 1.0);
    let mut rng = rand::thread_rng();
    let vals: Vec<_> = step.sample_iter(&mut rng).take(1000).collect();

    let mut group = c.benchmark_group("Aabb hit");
    //group.bench_function("Regular", |b| b.iter(|| bbox.hit_regular(&Ray::new(Vec3::new(0.0, 0.0, 0.0),Vec3::new(*vals.choose(&mut rng).unwrap(), *vals.choose(&mut rng).unwrap(), *vals.choose(&mut rng).unwrap())), &mut Interval::new(0.001,f64::INFINITY))));
    group.bench_function("Optimized", |b| {
        b.iter(|| {
            bbox.hit(
                &Ray::new(
                    Vec3::new(0.0, 0.0, 0.0),
                    Vec3::new(
                        *vals.choose(&mut rng).unwrap(),
                        *vals.choose(&mut rng).unwrap(),
                        *vals.choose(&mut rng).unwrap(),
                    ),
                ),
                &mut Interval::new(0.001, f64::INFINITY),
            )
        })
    });

    //group.bench_function("TM", | b | b.iter(|| t.hit_tm(&r, &mut i)));
    group.finish();

    let config: Config = Config { img_width: 300., img_height: 300., samples: 4, max_depth: 8, sky_color: Color::new(0.5, 0.5, 0.5), bvh_split_method: Some(BVHSplitMethod::Mid), };
    let json_camera: JSONCamera = JSONCamera { camera_center: Vec3::new(278., 278., -800.), camera_look_at: Vec3::new(278., 278., 0.), camera_fov_vertical: 40., camera_defocus_angle: 0., camera_focus_dist: 800. };
    
    let dragon = JSONElement::JSONObj(JSONObj{
        filepath: "input/obj/dragon.obj".to_string(),
        transpose: Some(Transpose {
            x: 220.,
            y: 0.,
            z: 200.
        }),
        rotate: Some(Rotate {
            theta_x: 0.,
            theta_y: -3.2,
            theta_z: 0.,
        }),
        scale: Some(Scale {
            x: 20.,
            y: 20.,
            z: 20.,
        }),
        material: raytracer::materials::Material::Lambertian(Lambertian::new(Color::new(0.7, 0.7, 0.7))),
    });

    let mut elements: Vec<JSONElement> = Vec::new();
    elements.push(dragon);

    let json_scene: JSONScene = JSONScene { 
        camera: json_camera, 
        elements,
    };

    let mut r = RenderIntegrator::new(json_scene, config);

    let mut group = c.benchmark_group("Render test");
    group.bench_function("Mid split", |b| {
        b.iter(|| {
            r.render()
        })
    });

    group.finish();

}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
