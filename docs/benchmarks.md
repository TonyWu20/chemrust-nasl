# Benchmark specification

Scope: benchmark work in the `chemrust-nasl` library crate.
Three parts, one example each in `chemrust-nasl/examples/`.
Shared code lives in `examples/bench_common/`.

## Part 1: grid-search completeness baseline (`grid_baseline`)

A brute-force voxel reference that the exact result is checked against.

- A voxel grid fills the Cartesian box of the model atoms.
- A grid point is a candidate when no atom is closer than r and at
  least one atom sits in the window [r, r + delta].
- Double candidates need two atoms in the window.
- Multi candidates need three or more.

A z-pinned variant of this scan runs inside `scalability` (Part 3):
the in-plane box is the surface area and the z range is the surface
+ vacuum band. This keeps the exact-versus-grid comparison on
surface area.

Completeness is measured against a fixed physical tolerance `delta`,
independent of the grid step. Per step the result records:

- scan time and candidate count (cost and redundancy),
- matched / missed sites (completeness),
- the largest site-to-candidate error (accuracy).

`--sweep N` runs the baseline at step, step/2, ..., step/2^(N-1) and
emits one block per step. This exposes the trade-off: a coarse step
misses sites, a fine step covers them but pays in candidate count and
scan time.

## Part 2: runtime benchmark with hardware record (`search_time`)

- Measure `search_sites` over repeated runs. Record wall time and the
  site count per class.
- Each run's site set is hashed. The `deterministic` field reports
  whether every run returned the identical site set.
- Write one JSON result file per run.

Every result file carries this context block:

- CPU model name, core count, and the rayon thread count used.
- Memory size.
- Operating system and kernel.
- Rust toolchain version, code git commit, and build profile.
- A note on frequency scaling when the host cannot pin clocks.

## Part 3: surface-area scalability run (`scalability`)

- Build an in-memory fcc slab with a fixed layer count and vacuum
  gap. The z axis is pinned to the surface + vacuum band, mirroring
  the app per-axis range flags. Only the surface-band atoms act as
  bond anchors.
- Scale the in-plane supercell, i.e. the surface area. This is the
  quantity that matters for high-throughput surface screening.
- Time `search_sites` on each in-plane size. Record atom count,
  anchor count, wall time, and the site count per class.
- `--grid` adds a z-pinned voxel-grid baseline at the same surface
  area, so the exact-versus-grid comparison is on surface area,
  not layer depth.
- The curve shows the exact engine scaling with surface area. The
  grid emits redundant candidates and misses sites at a coarse
  step, at the same surface area.

Slab model generation for any base uses `make_slab`, which chains
`chemrust-geometry` transforms: surface rotation, in-plane supercell,
layer replication, axis alignment, and a vacuum gap.
