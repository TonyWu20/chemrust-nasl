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
- Candidates cluster by 26-neighborhood.

Completeness is measured against a fixed physical tolerance `delta`,
independent of the grid step. Per step the result records:

- scan time and candidate count (cost side),
- matched / missed sites (completeness side),
- spurious clusters (redundancy side),
- the largest site-to-candidate error.

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

## Part 3: scalability run (`scalability`)

- Build an in-memory fcc slab family of growing atom count.
- Time `search_sites` on each family member.
- Record atom count, wall time, and the site count per class.
- The curve shows how the exact engine scales. Cost tracks the
  candidate-pair count (the local bond graph), so thick multi-layer
  slabs cost more per atom than a single surface.

Slab model generation for any base uses `make_slab`, which chains
`chemrust-geometry` transforms: surface rotation, in-plane supercell,
layer replication, axis alignment, and a vacuum gap.
