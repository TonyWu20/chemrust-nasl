//! Build vacuum slab models with `chemrust-geometry`.
//!
//! Pipeline: base structure -> SurfaceRotation(hkl) -> in-plane
//! supercell -> layer replication along c -> align axes -> vacuum
//! gap -> wrap fractional coordinates. The result writes as a
//! CASTEP `.cell` file, ready for the search and grid benchmarks.
//!
//! The (hkl) rotation suits cubic conventional cells.
//!
//! Usage:
//!   cargo run -p chemrust-nasl --example make_slab -- \
//!       fcc --species Pt --a 3.924 1 1 1 --layers 4 --super-cell 2 2 \
//!       --vacuum 12.0 --out Ptfcc111_4x22_vac12.cell
//!   cargo run -p chemrust-nasl --example make_slab -- \
//!       <base.cell|base.cif> 1 1 1 --layers 4 --super-cell 2 2 \
//!       --vacuum 12.0 --out slab.cell

mod bench_common;

use std::fs::write;
use std::path::Path;
use std::str::FromStr;

use castep_cell_fmt::{format::to_string_many_spaced, ToCellFile};
use castep_cell_io::{
    cell::{
        lattice_param::LatticeCart,
        positions::{PositionFracEntry, PositionsFrac},
        species::Species,
    },
    CellDocument, Lattice, Positions,
};
use chemrust_geometry::{slab::fcc_bulk, transform::SurfaceRotation, Structure};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mut base: Option<String> = None;
    let mut hkl: [i32; 3] = [0, 0, 0];
    let mut hkl_it = 0;
    let mut layers = 4;
    let mut super_cell: [usize; 2] = [2, 2];
    let mut sc_it = 0;
    let mut pending_super_cell = false;
    let mut vacuum = 12.0;
    let mut species = "Cu".to_string();
    let mut a_latt = 3.615;
    let mut out: Option<String> = None;
    let mut it = args.into_iter().skip(1);
    while let Some(a) = it.next() {
        match a.as_str() {
            "--layers" => {
                layers = it
                    .next()
                    .expect("--layers takes a value")
                    .parse()
                    .expect("number")
            }
            "--super-cell" => pending_super_cell = true,
            "--vacuum" => {
                vacuum = it
                    .next()
                    .expect("--vacuum takes a value")
                    .parse()
                    .expect("number")
            }
            "--species" => {
                pending_super_cell = false;
                species = it.next().expect("--species takes a value")
            }
            "--a" => {
                a_latt = it
                    .next()
                    .expect("--a takes a value")
                    .parse()
                    .expect("number")
            }
            "--out" => out = Some(it.next().expect("--out takes a value")),
            _ if base.is_none() => base = Some(a),
            _ if hkl_it < 3 => {
                hkl[hkl_it] = a.parse().expect("hkl index");
                hkl_it += 1;
            }
            _ if pending_super_cell && sc_it < 2 => {
                super_cell[sc_it] = a.parse().expect("supercell index");
                sc_it += 1;
                if sc_it == 2 {
                    pending_super_cell = false;
                }
            }
            _ => eprintln!("unexpected argument {a:?}"),
        }
    }
    if hkl_it != 3 || base.is_none() || out.is_none() {
        eprintln!("usage: <fcc|model> <h> <k> <l> [--layers N] [--super-cell NX NY] \
                   [--vacuum A] [--out PATH]");
        std::process::exit(1);
    }
    let base_arg = base.clone().expect("checked above");
    let out_path = out.clone().expect("checked above");

    let structure: Structure = if base_arg == "fcc" {
        let sym = castep_periodic_table::element::ElementSymbol::from_str(&species)
            .expect("a valid element symbol for --species");
        fcc_bulk(a_latt, sym)
    } else {
        bench_common::bench_model::load_structure(Path::new(&base_arg)).expect("load base model")
    };

    let slab = structure
        .transform(SurfaceRotation::new(hkl[0], hkl[1], hkl[2]))
        .apply()
        .supercell(super_cell[0], super_cell[1], 1)
        .replicate_along_c(layers)
        .align_axes()
        .add_vacuum_gap(vacuum)
        .wrap_frac_coords();

    let doc = to_cell_document(&slab);
    let text = to_string_many_spaced(&doc.to_cell_file());
    write(&out_path, text).expect("write slab file");
    let (la, lb, lc) = slab.require_cell().expect("lattice").lengths();
    println!(
        "wrote {out_path}: {} atoms, cell {la:.4} x {lb:.4} x {lc:.4} A",
        slab.num_atoms()
    );
}

fn to_cell_document(s: &Structure) -> CellDocument {
    let tensor = s
        .cell
        .as_ref()
        .expect("the slab carries a lattice")
        .tensor();
    let lattice = Lattice::Cart(LatticeCart {
        unit: None,
        a: tensor.column(0).into(),
        b: tensor.column(1).into(),
        c: tensor.column(2).into(),
    });
    let entries: Vec<PositionFracEntry> = s
        .species
        .iter()
        .zip(s.frac_coords.iter())
        .map(|(symbol, coord)| {
            PositionFracEntry::builder()
                .species(Species::Symbol(symbol.to_string()))
                .coord(coord.into_array())
                .build()
        })
        .collect();
    CellDocument::builder()
        .lattice(lattice)
        .positions(Positions::Frac(PositionsFrac { positions: entries }))
        .build()
        .expect("cell document builds")
}
