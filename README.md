![GitHub CI Status](https://img.shields.io/github/actions/workflow/status/scott223/raytracer/rust.yml?style=flat-square&logo=github)
[![dependency status](https://deps.rs/repo/github/scott223/raytracer/status.svg)](https://deps.rs/repo/github/scott223/raytracer)
[![lines count](https://img.shields.io/endpoint?url=https://ghloc.vercel.app/api/scott223/raytracer/badge?filter=.rs$,.toml$)](https://ghloc.vercel.app/scott223/raytracer?filter=.rs$,.toml$)


# Raytracer
Simple Raytracer for simple 3D scenes, made in Rust. Objective is to learn Rust, not make a very good raytracer.

## Ray tracing features
- monte carlo based, pixel by pixel raytracing with stochastic anti aliasing
- weigthed importance sampling
- various primitives (triangle, quad, sphere)
- various materials (Lambertian, Metal, Glass, ...)
- Bounded Hierarchy Volumes (BHV) tree to accelerate intersection tracing
- parallel multi-threaded rendering
- shadow acne prevention
- gamma correction
- configurable depth of focus
- scene and configuration input from JSON
- ...

## todo's
- jittering sampling

## Credits / references
- Raytracing in One Weekend (https://raytracing.github.io/books/RayTracingInOneWeekend.html)
- PBR book (https://www.pbr-book.org)
- Raytracer in Rust (https://github.com/dps/rust-raytracer)
- Raytracer in Rust (https://bheisler.github.io/post/writing-raytracer-in-rust-part-1/)

## License

Licensed under either of

 * Apache License, Version 2.0
   ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
 * MIT license
   ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

at your option.

## Contribution

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in the work by you, as defined in the Apache-2.0 license, shall be
dual licensed as above, without any additional terms or conditions.