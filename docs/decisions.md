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

- The `rhino` app-crate migration to the new stack is done. See
  the app-crate entries below. `elements_rs` still provides no
  CASTEP potential-name lookup. The seeding flow gets it from the
  `castep-periodic-table` lookup instead.
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

Note: the app-crate `chemrust-core` dependency named here was
removed in the next entry. The app now runs on `chemrust-geometry`
and the `castep-periodic-table` symbol table.

## App crate: off chemrust-core onto chemrust-geometry Structure

Date: 2026-07-10

`chemrust-core` is deprecated; `chemrust-geometry` is its successor. The
successor does not carry the `CrystalModel`/`CoreAtomData` trait pair.
Its model is the concrete `Structure` type (species, fractional
coordinates, `LatticeVectors`). The app crate's model layer now runs
on it:

- `execution/model.rs` converts between `Structure` and the two file
  formats. `cell_to_structure` reads a `castep-cell-io` 0.7
  `CellDocument` (lattice plus per-entry species).
  `cif_to_structure` reads a `crystal-cif-io` 0.2 `CoreCif`. Both
  export directions build CIF text and CASTEP text from a `Structure`
  plus one new adsorbate site.
- `SearchJob`, `ExportFormat`, and `RhinoExport` all take a
  `Structure`. The `.cell` I/O moved from `castep-cell-io` 0.2.18 to
  0.7.0 with `castep-cell-fmt` 0.3.0 parsing and writing.
- The seeding seam (`castep-seeding` traits, `RootJobs`,
  `build_all`) sat on the 0.2.18 stack through a renamed dependency
  `castep_cell_io_v02` at this point. The next entry records its
  removal.

Verification: the round-trip test searches `SAC_GDY_V.cell` through
the 0.7 parser (24 single, 51 double, 19 multi sites at 2.5 A),
exports one site as CIF, re-parses it, and re-runs the search. The
counts match. The lib crate's 73/93/26 baseline is unchanged.

## castep-seeding 0.2.0: the seeding seam moves to the 0.7 stack

Date: 2026-07-10

The `castep-seeding` repo was migrated in its own `dev` branch
(commit `1264694`, pushed): published `castep-cell-io 0.7.0` plus
`castep-cell-fmt 0.3.0` for parsing and writing, the param model
inlined from `castep-param-io` into `src/param`, a local
`EnergyCutoff` module replacing the 0.2.18 one, and `CellDocument`
helpers (`cell_elements`, `total_spin`, `species_params_for`,
`cell_to_string`) carrying over the 0.2.18 semantics. Its
`CastepTask` gained a `SinglePoint` variant under
`#[non_exhaustive]`.

In this repo the app crate now consumes it:

- `castep-seeding = { git = "https://github.com/TonyWu20/castep-seeding",
  rev = "1264694..." }` pins the pushed commit. The dependency is
  still a git dep because no 0.7-stack release is published. It
  swaps to a registry dep when one lands.
- The `castep_cell_io_v02` alias and the `castep-cell-io 0.2.18`
  and `castep-param-io` dependencies are removed.
  `CastepTask`, `EnergyCutoff`, `CastepParam`, and
  `parse_cell_doc_from_path` all come from `castep-seeding` now.
- `ExportFile` drops its `T: Display` bound. The 0.7 `CellDocument`
  has no `Display` impl, and the seeding flow serializes it through
  `castep_seeding::cell_to_string`.
- Task matches gain a wildcard arm because `CastepTask` is
  `#[non_exhaustive]`.

Verification: `cargo test -p chemrust-nasl` passes 6 tests
(including the 73/93/26 search baseline) and `cargo test -p rhino`
passes 4 (including the 24/51/19 CIF round-trip). `Cargo.lock`
carries no `chemrust-core` and no `castep-cell-io 0.2.18` entry
anymore.
