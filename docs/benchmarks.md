# Benchmark specification

Scope: benchmark work in the `chemrust-nasl` library crate.
Two parts, one example each in `chemrust-nasl/examples/`.

## Part 1: grid-search completeness baseline

A brute-force reference that checks the exact search results.

- A voxel grid fills the Cartesian search box.
- The step size is configurable. The default is 0.05 A.
- A grid point is a candidate when no atom is closer than r and at
  least one atom sits in the window [r, r + delta]. The default
  delta is 0.1 A.
- Double candidates need two atoms in the window.
- Multi candidates need three or more atoms in the window.
- Metrics against one `SearchConfig` run:
  - match rate: the share of exact sites with a candidate within one
    grid step.
  - miss rate: candidate clusters farther than one grid step from any
    exact site.
  - the largest distance error between a site and its nearest
    candidate.

## Part 2: runtime benchmark with hardware record

- Measure `search_sites` on structures of growing atom count.
- Record wall time and the site count per class.
- Write one JSON result file per run.

Every result file carries this context block:

- CPU model name, core count, and the rayon thread count used.
- Memory size.
- Operating system and kernel.
- Rust toolchain version, code git commit, and build profile.
- A note on frequency scaling when the host cannot pin clocks.

Implementation home: one example per part in `chemrust-nasl/examples/`.
