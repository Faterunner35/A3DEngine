# A3DEngine

An experimental C++20 3D engine built from the mathematics up. The core is a small, dependency-free math library, a data-oriented ECS for transforms, and geometric algebra (PGA, Cl(3,0,1)) for points, lines and planes.

The project is a research platform. Each subsystem is written from first principles, checked against a reference implementation, and measured.

---

## What is here now

### Math core (`SegmentTest/`)
| File | What it does |
|---|---|
| `Geomerty.h` | `Vector<N,T>`: dot, cross, robust length, Gram–Schmidt orthonormalisation |
| `Matrix3x3.hpp`, `Matrxi4x4.hpp` | 3×3 / 4×4 matrices, determinant, adjugate inverse |
| `Affine4X3.hpp` | Affine transform stored as 3×4 (the `0 0 0 1` row is implied): composition in 36 multiplies, inverse via cofactors with a relative (Hadamard) singularity test, rigid inverse, normal matrix `sign(det)·cof(A)` |
| `Quaternion.hpp`, `QuatSquad.hpp`, `QuatInterpolationPrecomputed.hpp` | Quaternions, slerp/squad interpolation |
| `TransfromMatrix.hpp` | `MakeTransform(T, R, S)` and its exact inverse built from TRS without inverting a matrix |
| `EulerAngles.hpp`, `PrimitiveTransform.hpp` | Euler angles, elementary rotations |

### Geometry
`Plane.h`, `Ray.h`, `Sphere.h`, `AABB.h`, `GeomIntersector.h`, `ManifoldMesh.h`: primitives, intersections, and a manifold triangle mesh with edge/triangle adjacency.

### Geometric algebra
`WedgeProduct.h` is a `Multivector` for projective geometric algebra Cl(3,0,1) with wedge `^`, regressive product `&` (meet/join), geometric product `*` and duality. A plane through three points is just `R ^ S ^ T`. A line where two planes meet is `P & Q`.

### ECS: entities and transforms
`Entity.hpp` and `TranformManager.hpp`:

- **`EntityPool`**: a fixed-capacity pool of `{index, generation}` handles. `Create` and `Destroy` are O(1). Freed slots are reused in FIFO order, so a stale handle is always detected.
- **`TransformComponent`**: dense structure-of-arrays storage behind a sparse set.
  - **Input:** local TRS (position, quaternion, scale) relative to the parent.
  - **Output:** world transform as an exact affine 3×4 matrix `(A, t)`, computed as `W = W_parent · T·R·S`.
  - **Hierarchy:** arrays are kept in topological order (sorted by depth with a counting sort), so `Update()` is **one linear pass** with no recursion and no pointer chasing.
  - **No world TRS is stored.** With non-uniform scale, a parent's scale combined with a child's rotation produces shear, which TRS cannot represent (see test 1.2: the axes end up 126.87° apart). The world matrix is kept exact. Derived data (inverse, normal matrix, polar decomposition, bounding radius) is computed by whoever consumes it.

---

## Build

Requires CMake ≥ 3.20 and a C++20 compiler (MSVC 2022, GCC ≥ 11, Clang ≥ 14).

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
ctest --test-dir build -C Release          # correctness tests only
./build/engine_tests                       # tests + benchmark + hierarchy stress
./build/engine_tests --bench               # tests + benchmark only
./build/engine_tests --stress              # tests + hierarchy stress only
#   Windows: build\Release\engine_tests.exe
```

### Targets

| Target | Source | What it is |
|---|---|---|
| `a3d_core` | everything in `SegmentTest/` except entry points | static library: the engine |
| `engine_tests` | `SegmentTest/main.cpp` | correctness tests, benchmark, hierarchy stress for `EntityPool` + `TransformComponent` |
| `pga_demo` | `SegmentTest/PgaDemo.cpp` | PGA Cl(3,0,1) demo: plane through three points, geometric product, meet |

CMake is the only build description. There is no hand-written `.sln`.

### Visual Studio 2022

1. *File → Open → Folder…* → the repository root. VS reads `CMakeLists.txt` and `CMakePresets.json`.
2. Pick `x64 Debug` or `x64 Release` in the configuration box, `engine_tests.exe` in the startup item box, F5.

**Adding a source file.** In Solution Explorer switch to *Folder View* (the "Switch between solutions and available views" button), right-click `SegmentTest` → *Add → New File…* (or *Add → Existing Item…*). Nothing has to be edited in `CMakeLists.txt`: sources are collected with `file(GLOB_RECURSE … CONFIGURE_DEPENDS)`, so on the next build CMake notices the new file, regenerates and compiles it into `a3d_core`. Subfolders of `SegmentTest/` work too. Deleting or renaming a file is picked up the same way.

A new file with its own `main()` is a new program, not part of the engine: add it to `A3D_ENTRY_POINTS` and give it an `add_executable` in `CMakeLists.txt`.

If you prefer a classic solution (for example for the profiler): `cmake --preset vs2022` creates `out/build/vs2022/A3DEngine.sln`. That solution is generated output, so don't add files through its *Project → Add New Item* (it would put them under `out/`); create them in `SegmentTest/` and the `ZERO_CHECK` project picks them up on the next build.

| Option | Effect |
|---|---|
| `-DA3D_NATIVE=ON` | `-march=native` / `/arch:AVX2` |
| `-DA3D_SANITIZERS=ON` | AddressSanitizer + UBSan (GCC/Clang) |

---

## Benchmark: `TransformComponent::Update()`

Measured on a cloud VM (2 vCPU Intel Xeon @ 2.1 GHz, GCC 13.3, `-O3`), so treat the numbers as relative rather than absolute. Every frame recomputes **every** world matrix. Scenes are built in shuffled order, so the first frame must also rebuild the hierarchy order.

| N | hierarchy | first frame (with reorder) | per frame | ns / transform | share of a 60 FPS frame |
|---:|---|---:|---:|---:|---:|
| 10 000 | flat | 0.16 ms | **0.07 ms** | 6.9 | 0.4 % |
| 10 000 | random forest | 1.24 ms | **0.13 ms** | 12.8 | 0.8 % |
| 10 000 | skeletons (chains of 32) | 0.53 ms | **0.14 ms** | 13.7 | 0.8 % |
| 50 000 | random forest | 5.95 ms | **0.82 ms** | 16.4 | 4.9 % |
| 100 000 | flat | 1.51 ms | **0.84 ms** | 8.4 | 5.0 % |
| 100 000 | random forest | 10.6 ms | **2.54 ms** | 25.4 | 15 % |
| 100 000 | skeletons | 9.8 ms | **1.71 ms** | 17.1 | 10 % |
| 1 000 000 | flat | 17.5 ms | 17.4 ms | 17.4 | 104 % |
| 1 000 000 | random forest | 707 ms | 55.7 ms | 55.7 | 334 % |

Compared with a classic pointer-based node tree (one heap allocation per node, children by pointer, recursive update) on identical data, with bit-identical results:

| N | TransformComponent | pointer tree | speed-up |
|---:|---:|---:|---:|
| 10 000 | 0.13 ms | 0.93 ms | ×7.3 |
| 100 000 | 3.0 ms | 17.0 ms | ×5.7 |
| 1 000 000 | 67.6 ms | 343.6 ms | ×5.1 |

Memory per transform is 200 bytes: 80 for local TRS, 96 for the world 3×4 matrix, and 24 for handles and indices.

**Reading the benchmark numbers.**
- **Up to ~100 k transforms:** a full update costs 1–15 % of a 60 FPS frame on a single core.
- **Beyond that:** the working set (200 B × N) no longer fits in cache and the update becomes memory-bound. Random parents cost more than flat scenes because the parent's matrix is a non-sequential read.
- **Reorder cost:** the order rebuild is paid only when the hierarchy changes, not every frame.
- **Next steps:** dirty-subtree updates and splitting `Update()` by depth level across threads.

### Hierarchy stress (`--stress`)

| Scenario | Result |
|---|---|
| Chain 20 000 levels deep, each node rotated 2π/1000 (exact answer: a regular 1000-gon) | `Update` 0.6 ms; error 2.6e-13 in rotation, 4.1e-11 in translation |
| One root with 100 000 children, then the root is destroyed | `Update` 2.6 ms; destroy + GC + `Update` 7.2 ms; all 100 000 children become roots |
| Live world, 50 000 entities, every frame: 0.5 % die, 0.5 % spawn, 0.5 % reparent, 10 % animate | frame p50 4.2 ms, p99 6.5 ms; bit-identical to the reference |
| Heavy mode, 50 000 entities: 2 % die, 2 % spawn, 5 % reparent, 30 % animate | frame p50 6.7 ms, p99 11.3 ms |
| 100 000 entities: one `SetParent` that breaks the order | next `Update` 8.1 ms instead of 2.9 ms (full reorder) |
| Live world, **200 000** entities, 0.5 % die / spawn / reparent, 10 % animate | frame p50 49 ms, p99 80 ms; bit-identical to the reference |
| **200 000** entities, animation only (no structural change) | frame p50 31 ms (`Update` 25 ms) |
| Live world, **500 000** entities, 0.5 % die / spawn / reparent, 10 % animate | frame p50 156 ms, p99 239 ms; bit-identical to the reference |
| **500 000** entities, heavy mode (2 % / 2 % / 5 % / 30 %) | frame p50 288 ms; 1.48 M reparents over 60 frames |

Known costs these tests expose:
- `SetParent` walks up the ancestors to reject cycles, so it is O(depth). Building a 20 000-deep chain node by node is quadratic (1.3 s).
- Any `Destroy`, and any reparent onto a node further right in the array, invalidates the whole order. In a world with deaths every frame, the reorder is paid every frame.

- Beyond ~100 000 entities a full `Update` becomes memory-bound. When the hierarchy is built incrementally, the arrays stay in insertion order rather than depth order, so reading the parent's matrix is a random access across 40–95 MB. At 200 000 entities this alone costs about 25 ms per frame, even when only 10 % of the transforms changed.

Planned fixes: keep a depth per entity for O(1) cycle rejection in the common case; make the reorder incremental (move the affected subtree only); recompute only dirty subtrees; order siblings so that children follow their parent closely in memory; update level by level across threads.

---

## Roadmap

### Progressive geometry
Continuous level of detail on top of `ManifoldMesh`:
- progressive meshes (Hoppe): *edge collapse / vertex split* sequences with quadric error metrics;
- streaming refinement: a coarse base mesh plus a stream of vertex splits;
- view-dependent refinement driven by the world transform (screen-space error from the spectral norm of `A`);
- keeping the mesh manifold and preserving attributes (normals, UVs) across collapses.

### Grassmann (exterior) algebra
Growing `WedgeProduct.h` into a full exterior-algebra layer:
- a graded `k`-vector API (points, lines, planes as 1-, 2-, 3-vectors) with wedge, regressive product and Hodge duality;
- Plücker coordinates for lines, and meet/join for intersections without special cases;
- determinants, areas and volumes as wedge products, with orientation tests for robust geometry;
- motors (PGA rotors and translators) as an alternative transform representation, compared against the affine 3×4 path in the same benchmark.

### Engine
- polar decomposition `A = R·P` for world-space interpolation (motion blur, physics, editor gizmos);
- a geometry component (local-space meshes, world bounds via the spectral norm);
- a SIMD 3×4 composition path and a multithreaded level-by-level update.

---

## Layout

```
A3DEngine/
├── CMakeLists.txt
├── SegmentTest.sln            # original Visual Studio solution
└── SegmentTest/
    ├── main.cpp               # ECS tests, benchmark, stress (engine_tests)
    ├── SegmentTest.cpp        # PGA demo (not built yet)
    ├── Entity.hpp             # EntityPool
    ├── TranformManager.hpp    # TransformComponent
    ├── Affine4X3.hpp, TransfromMatrix.hpp, Quaternion.hpp, Matrix3x3.hpp, ...
    ├── WedgeProduct.h         # PGA multivector
    └── ManifoldMesh.h, Plane.h, Ray.h, Sphere.h, AABB.h, ...
```
