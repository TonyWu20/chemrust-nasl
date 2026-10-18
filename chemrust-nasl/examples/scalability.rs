//! Scalability run: time the exact search on an in-memory fcc slab
//! family of growing atom count. This curve supports the complexity
//! claim (kd-tree build O(n log n), per-atom queries O(log n + k),
//! O(1) per closed-form intersection).
//!
//! Each family member is a cubic (hkl) surface: an fcc bulk primitive
//! rotated to the surface, supercell in-plane, layers replicated
//! along c, and a vacuum gap. The result is JSON plus a table.
//!
//! Usage:
//!   cargo run -p chemrust-nasl --release --example scalability -- \
//!       [--species Pt] [--a 3.924] [--hkl 1 1 1] [--vacuum 12.0] \
//!       [--bondlength 2.5] [--runs 3] [--out FILE.json]

mod bench_common;

use std::fs::write;
use std::path::PathBuf;
use std::str::FromStr;
use std::time::Instant;

use castep_periodic_table::element::ElementSymbol;
use chemrust_geometry::{slab::fcc_bulk, transform::SurfaceRotation, Structure};
use chemrust_nasl::SearchConfig;
use nalgebra::Point3;
use serde::Serialize;

/// Default family: (layers, in-plane x, in-plane y).
const SERIES: &[(usize, usize, usize)] = &[
    (4, 1, 1),
    (4, 2, 2),
    (4, 3, 3),
    (4, 4, 4),
    (8, 2, 2),
    (8, 3, 3),
    (8, 4, 4),
];

#[derive(Serialize)]
struct Row {
    layers: usize,
    in_plane: String,
    atom_count: usize,
    runs_ms: Vec<f64>,
    multi: usize,
    single: usize,
    double: usize,
}

#[derive(Serialize)]
struct ResultJson {
    benchmark: &'static str,
    species: String,
    hkl: [i32; 3],
    bondlength: f64,
    vacuum: f64,
    rows: Vec<Row>,
    context: bench_common::hardware_spec::Context,
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mut species = "Pt".to_string();
    let mut a_latt = 3.924_f64;
    let mut hkl: [i32; 3] = [1, 1, 1];
    let mut hkl_it = 0;
    let mut vacuum = 12.0_f64;
    let mut bondlength = 2.5_f64;
    let mut runs = 3usize;
    let mut out: Option<PathBuf> = None;
    let mut it = args.into_iter().skip(1);
    while let Some(a) = it.next() {
        match a.as_str() {
            "--species" => species = it.next().expect("--species takes a value"),
            "--a" => a_latt = it.next().expect("--a takes a value").parse().expect("number"),
            "--vacuum" => vacuum = it.next().expect("--vacuum takes a value").parse().expect("number"),
            "--bondlength" => bondlength = it.next().expect("--bondlength takes a value").parse().expect("number"),
            "--runs" => runs = it.next().expect("--runs takes a value").parse().expect("number"),
            "--out" => out = Some(PathBuf::from(it.next().expect("--out takes a value"))),
            _ if hkl_it < 3 => {
                hkl[hkl_it] = a.parse().expect("hkl index");
                hkl_it += 1;
            }
            _ => eprintln!("unexpected argument {a:?}"),
        }
    }
    if hkl_it != 3 {
        eprintln!("usage: [--species S] [--a A] <h> <k> <l> [--vacuum V] \
                   [--bondlength R] [--runs N] [--out FILE]");
        std::process::exit(1);
    }

    let sym = ElementSymbol::from_str(&species).expect("a valid element symbol for --species");
    let mut rows = Vec::with_capacity(SERIES.len());
    for &(layers, nx, ny) in SERIES {
        let structure = slab(layers, nx, ny, &sym, a_latt, hkl, vacuum);
        let points = bench_common::bench_model::cart_points(&structure);
        let to_check: Vec<(usize, Point3<f64>)> =
            points.iter().enumerate().map(|(i, p)| (i, *p)).collect();
        let config = SearchConfig::new(&to_check, &points, bondlength);
        let report = config.search_sites();
        let mut runs_ms = Vec::with_capacity(runs);
        for _ in 0..runs {
            let t = Instant::now();
            config.search_sites();
            runs_ms.push(t.elapsed().as_secs_f64() * 1e3);
        }
        let multi = report.points().map(|v| v.len()).unwrap_or(0);
        let single = report.viable_single_points().map(|v| v.len()).unwrap_or(0);
        let double = report.viable_double_points().map(|v| v.len()).unwrap_or(0);
        let mean_ms = runs_ms.iter().sum::<f64>() / runs.max(1) as f64;
        rows.push(Row {
            layers,
            in_plane: format!("{nx}x{ny}"),
            atom_count: points.len(),
            runs_ms,
            multi,
            single,
            double,
        });
        println!(
            "{species}{hkl:?} layers={layers} in-plane={nx}x{ny} atoms={:3}  \
             multi/single/double={}/{}/{}  ~{:.1} ms",
            points.len(),
            multi,
            single,
            double,
            mean_ms
        );
    }

    let json = ResultJson {
        benchmark: "scalability",
        species: species.clone(),
        hkl,
        bondlength,
        vacuum,
        rows,
        context: bench_common::hardware_spec::Context::capture(rayon::current_num_threads()),
    };
    let text = serde_json::to_string_pretty(&json).expect("serialize result");
    match out {
        Some(path) => write(path, text).expect("write result file"),
        None => println!("{text}"),
    }
}

/// Build one slab in memory: fcc bulk -> (hkl) surface -> in-plane
/// supercell -> layer replication -> align axes -> vacuum -> wrap.
fn slab(
    layers: usize,
    nx: usize,
    ny: usize,
    sym: &ElementSymbol,
    a_latt: f64,
    hkl: [i32; 3],
    vacuum: f64,
) -> Structure {
    fcc_bulk(a_latt, *sym)
        .transform(SurfaceRotation::new(hkl[0], hkl[1], hkl[2]))
        .apply()
        .supercell(nx, ny, 1)
        .replicate_along_c(layers)
        .align_axes()
        .add_vacuum_gap(vacuum)
        .wrap_frac_coords()
}
