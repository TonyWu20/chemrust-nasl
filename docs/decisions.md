# Decision records

Session of 2026-10-12. This document records the decisions taken when
migrating the workspace to the successor crate stack.

## Confirmed decisions

- **Migration scope.** The deprecated `chemrust-core` dependency is
  replaced by `chemrust-geometry` only in the library crate
  `chemrust-nasl`. The binary crate `rhino` (chemrust-nasl-app)
  keeps the old stack until its I/O dependencies move.
- **`elements_rs` features.** Declared with
  `default-features = false`, matching `crystal-cif-io`.
- **Test dev-deps.** The library tests use published registry
  `castep-cell-io 0.7.0` and `castep-cell-fmt 0.3.0`. Local path
  deps on sibling repos do not survive a fresh clone or CI.
- **Repo hygiene.** `scanner_test_models/` and `Cargo.lock` are
  tracked. The runtime seed output directories `V_1.5_CuTP001/`
  and `V_1.5_Ag_111_v20_1_PFGEV_ot/` are ignored.

## Behavior notes

- `Visualize::element_by_cn_number` keeps the legacy
  `castep-periodic-table` mapping: CN 0 maps to H, CN 1 falls
  back to W (the castep enum had no discriminant 1), CN 2 to 103
  maps to the element with that atomic number, and CN 104 or
  more maps to Np. See `visualize.rs` for the full rationale.
- `Atom.coord` is a Cartesian `Point3<f64>` in Angstrom. The
  legacy type was `chemrust_core::CoordData`
  (Fractional/Cartesian enum). All callers in this repo pass
  Cartesian coordinates.

## Open items

- Migrate `rhino` to the new stack: `castep-cell-io 0.7.0`,
  `crystal-cif-io dev-0.2` (elements_rs based), and a successor
  to `castep-seeding`. This rewrites `format_loader`, `export`,
  `seeding`, and the element plumbing. `elements_rs` provides no
  CASTEP potential-name lookup, so the seeding flow needs that
  data from another source.
- Refresh the README: the binary is `rhino`, the config example
  lacks `build_seed` and `potential_dir`, and several claims
  (browser use, binary download site) do not match the repo.
- Decide whether to keep the two `target.<triple>.linker` keys
  in the app manifest. Cargo 1.99 warns that it ignores them.
- `RunMode` is collected in the interactive UI, but no code
  consumes it.
- `ask_potential_dir` offers `e.g.: /path/Potentials` as its
  default. Accepting it writes the prefix into the path.
- Add a CI workflow.
