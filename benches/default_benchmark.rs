use criterion::{criterion_group, criterion_main, Criterion};
use rand::{rngs::SmallRng, SeedableRng};
use raytracer::vec3::Vec3;
use raytracer::elements::{Triangle, Hittable};
use raytracer::materials::{Lambertian, Material};
use raytracer::color::Color;
use raytracer::ray::Ray;
use raytracer::interval::Interval;

pub fn criterion_benchmark(c: &mut Criterion) {

    // random unit sphere
    let mut rng = SmallRng::seed_from_u64(223);
    c.bench_function("random_unit_sphere", |b| b.iter(|| Vec3::new_random_unit_sphere(&mut rng)));

    // triangle hit
    let m1: Material = Material::Lambertian(Lambertian::new(Color::new(1.0, 1.0, 1.0)));
    let t: Triangle = Triangle::new(
        Vec3::new(-2.0, -2.0, -5.0),
        Vec3::new(-2.0, 2.0, -5.0),
        Vec3::new(2.0, 2.0, -5.0),
        m1,
    );

    let r: Ray = Ray::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, -1.0));
    let mut i = Interval::new(0.001, f64::MAX);

    c.bench_function("triangle_hit", | b | b.iter(|| t.hit(&r, &mut i)));

}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);