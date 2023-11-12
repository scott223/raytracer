use criterion::{criterion_group, criterion_main, Criterion};
use rand::seq::SliceRandom;
use rand_distr::Distribution;
use rand_distr::Uniform;
use raytracer::bvh::aabb::Aabb;

use raytracer::interval::Interval;
use raytracer::render::ray::Ray;
use raytracer::linalg::vec3::Vec3;

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
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
