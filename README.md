![GitHub CI Status](https://img.shields.io/github/actions/workflow/status/scott223/raytracer/rust.yml?style=flat-square&logo=github)
[![dependency status](https://deps.rs/repo/github/scott223/raytracer/status.svg)](https://deps.rs/repo/github/scott223/raytracer)
[![lines count](https://img.shields.io/endpoint?url=https://ghloc.vercel.app/api/scott223/raytracer/badge?filter=.rs$,.toml$)](https://ghloc.vercel.app/scott223/raytracer?filter=.rs$,.toml$)


# Raytracer
Simple Raytracer for simple 3D scenes, made in Rust. Objective is to learn Rust, not make a very good raytracer.

## Ray tracing features
- Monte Carlo raytracer with adaptive sampling (based on z-test confidence interval) 
- Anti aliasing using Sobol sequence and Mitchell Netravali filter
- Weigthed importance sampling using Probability Distribution Functions and attractors (e.g. lights)
- Various primitives (triangle, quad, sphere)
- Various materials (Lambertian, Metal, Glass, ...)
- Bounded Hierarchy Volumes (BHV) tree, constructed using Surface Area Heuristic (SAH), to accelerate intersection tracing
- Parallel multi-threaded rendering using Rayon
- Gamma correction
- Configurable depth of focus
- Scene and configuration input from JSON
- ...

![Render example](https://github.com/scott223/raytracer/blob/main/render1200.png?raw=true)
*Render of different objects (cube, sphere, .obj file import (100k triangles) and materials (lambertian, glass, fuzzy metal). This render of 1200x1200 pixels with 2000 samples per pixel took about 30 minutes to complete on MacBook Air M2*

## todo's
- BRDF materials

## Rust learnings applied
- modules, imports, external crates
- ownership/references/borrow, explicit lifetimes
- methods, enums, structs
- generic methods
- traits
- basic iterators, match statement
- Result, Error, Option
- reading files, saving files
- logging, .env
- tests, benchmarking, flamegraph
- VSCode debugger, breakpoints and memory inspector
- Git(hub)

## Credits / references
- Raytracing in One Weekend (https://raytracing.github.io/books/RayTracingInOneWeekend.html)
- PBR book (https://www.pbr-book.org)
- Raytracer in Rust (https://github.com/dps/rust-raytracer)
- Raytracer in Rust (https://bheisler.github.io/post/writing-raytracer-in-rust-part-1/)

## License

Licensed under MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

## Contribution

Unless you explicitly state otherwise, any contribution intentionally submitted for inclusion in the work by you shall be licensed as above, without any additional terms or conditions.