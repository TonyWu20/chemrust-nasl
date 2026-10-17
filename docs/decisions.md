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

## App crate: crystal-cif-io 0.2.0 and the CIF model bridge

Date: 2026-07-10

The app crate (`chemrust-nasl-app`, binary `rhino`) now builds against the
published `crystal-cif-io 0.2.0` instead of the git dev 0.1.2. The 0.2.0
line split the generic `DataBlock` world into a typed `CoreCif` layer:
`parse_core_cif` returns `CoreCifDocument { blocks: Vec<CoreCif> }`, and
the 0.1.x `chemrust` bridge feature is gone.

The search engine still consumes `chemrust-core` traits
(`CrystalModel`, `CoreAtomData`). The new `execution/cif_model.rs` module
is the bridge between the two:

- `CifModel` owns a `LatticeVectors` and a `CifAtomData` (castep
  `ElementSymbol` symbols, `CoordData` coordinates, labels) read from a
  `CoreCif` cell items plus `atom_site` rows.
- CIF cell items carry Angstrom lengths and degree angles. Missing
  angles default to 90. Rows need fractional or Cartesian coordinates
  and a castep-known symbol.
- `CifExport` (in `execution/export.rs`) is the export item for the CIF
  path. It builds one `CoreCif` data block from the model atoms plus
  the new adsorbate site and writes it with `write_core_cif`. The rows
  carry labels, symbols, and coordinates only.

Verification: the round-trip test in `execution/export.rs` searches
`SAC_GDY_V.cell` (24 single, 51 double, 19 multi sites at 2.5 A with
z in [0.4, 0.5]), exports one site as CIF text, re-parses it, and
re-runs the search. The counts match the cell-model search exactly.

Open item: the lib crate (`chemrust-nasl`) is migrated to
`chemrust-geometry` and `elements_rs` (Decision A above). The app
crate still depends on `chemrust-core` and `castep-periodic-table`
through the search engine and this bridge. A full de-castep of the
search engine is a separate, larger task.
