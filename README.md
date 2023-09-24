# Raytracer
Simple Raytracer for simple 3D scenes, made in Rust. Objective is to learn Rust, not make a very good raytracer.

## Ray tracing features
- pixel by pixel raytracing with stochastic anti aliasing
- various primitives (sphere, quadrilateral)
- various materials (Lambertian, Metal, Glass, ...)
- Bounded Volume Hierarchy (BHV) tree to accelerate intersection tracing
- parallel multi-threaded rendering
- shadow acne prevention
- gamma correction
- configurable depth of focus
- scene and configuration input from JSON
- ...

## Credits / references
- Raytracing in One Weekend (https://raytracing.github.io/books/RayTracingInOneWeekend.html)
- Raytracer in Rust (https://github.com/dps/rust-raytracer)
- Raytracer in Rust (https://bheisler.github.io/post/writing-raytracer-in-rust-part-1/)