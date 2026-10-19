//! Runtime benchmark and determinism check, part 2 of
//! `docs/benchmarks.md`.
//!
//! Measures `SearchConfig::search_sites` across repeated runs. Each
//! run's site set is hashed. The `deterministic` field reports
//! whether every run returned the identical site set. This is the
//! evidence for the determinism claim.
//!
//! Usage:
//!   cargo run -p chemrust-nasl --release --example search_time -- \
//!       <model.cell|model.cif> <bondlength> [--runs N] [--out FILE.json]

mod bench_common;

use std::cmp::Ordering;
use std::fs::write;
use std::path::PathBuf;
use std::time::Instant;

use chemrust_nasl::{SearchConfig, SearchReports};
use nalgebra::Point3;
use serde::Serialize;

const DEFAULT_RUNS: usize = 5;

#[derive(Serialize)]
struct SiteCounts {
    multi: usize,
    single: usize,
    double: usize,
}

#[derive(Serialize)]
struct ResultJson {
    benchmark: &'static str,
    model: String,
    bondlength: f64,
    atom_count: usize,
    runs_ms: Vec<f64>,
    site_counts: SiteCounts,
    /// FNV-1a over the sorted site coordinates of each run.
    run_hashes: Vec<String>,
    /// True when every run returned the same site set.
    deterministic: bool,
    context: bench_common::hardware_spec::Context,
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mut model = None;
    let mut bondlength: Option<f64> = None;
    let mut runs = DEFAULT_RUNS;
    let mut out: Option<PathBuf> = None;
    let mut it = args.into_iter().skip(1);
    while let Some(a) = it.next() {
        match a.as_str() {
            "--runs" => runs = it.next().expect("--runs takes a value").parse().expect("number"),
            "--out" => out = Some(PathBuf::from(it.next().expect("--out takes a value"))),
            _ if model.is_none() => model = Some(PathBuf::from(a)),
            _ if bondlength.is_none() => {
                bondlength = Some(a.parse().expect("bondlength in Angstrom"))
            }
            _ => eprintln!("unexpected argument {a:?}"),
        }
    }
    let model_path = model.expect("usage: <model> <bondlength> [--runs N] [--out FILE]");
    let r = bondlength.expect("usage: <model> <bondlength> [--runs N] [--out FILE]");

    let structure = bench_common::bench_model::load_structure(&model_path).expect("load model");
    let points = bench_common::bench_model::cart_points(&structure);
    let to_check: Vec<(usize, Point3<f64>)> =
        points.iter().enumerate().map(|(i, p)| (i, *p)).collect();
    let config = SearchConfig::new(&to_check, &points, r);

    // Warm-up run carries the site counts (deterministic input).
    let report = config.search_sites();
    let counts = site_counts(&report);
    let mut runs_ms = Vec::with_capacity(runs);
    let mut run_hashes = Vec::with_capacity(runs);
    for _ in 0..runs {
        let t = Instant::now();
        let rep = config.search_sites();
        runs_ms.push(t.elapsed().as_secs_f64() * 1e3);
        run_hashes.push(site_set_hash(&rep));
    }
    let deterministic = !run_hashes.is_empty() && run_hashes.iter().all(|h| h == &run_hashes[0]);

    let json = ResultJson {
        benchmark: "search_time",
        model: model_path.to_string_lossy().to_string(),
        bondlength: r,
        atom_count: points.len(),
        runs_ms,
        site_counts: counts,
        run_hashes,
        deterministic,
        context: bench_common::hardware_spec::Context::capture(rayon::current_num_threads()),
    };
    let text = serde_json::to_string_pretty(&json).expect("serialize result");
    match out {
        Some(path) => write(path, text).expect("write result file"),
        None => println!("{text}"),
    }
}

fn site_counts(report: &SearchReports) -> SiteCounts {
    SiteCounts {
        multi: report.points().map(|v| v.len()).unwrap_or(0),
        single: report.viable_single_points().map(|v| v.len()).unwrap_or(0),
        double: report.viable_double_points().map(|v| v.len()).unwrap_or(0),
    }
}

/// FNV-1a over the sorted site coordinates of a report.
fn site_set_hash(report: &SearchReports) -> String {
    let bytes = site_set_bytes(report);
    let mut hash: u32 = 0x811c_9dc5;
    for b in bytes {
        hash ^= u32::from(b);
        hash = hash.wrapping_mul(0x0100_0193);
    }
    format!("{hash:08x}")
}

/// Little-endian bytes of every site coordinate, in sort order.
fn site_set_bytes(report: &SearchReports) -> Vec<u8> {
    let mut pts: Vec<Point3<f64>> = Vec::new();
    if let Some(v) = report.points() {
        pts.extend(v.iter().map(|p| p.point()));
    }
    if let Some(v) = report.viable_single_points() {
        pts.extend(v.iter().map(|p| p.point()));
    }
    if let Some(v) = report.viable_double_points() {
        pts.extend(v.iter().map(|p| p.point()));
    }
    pts.sort_by(|a, b| sort_key(a).partial_cmp(&sort_key(b)).unwrap_or(Ordering::Equal));
    let mut out = Vec::with_capacity(pts.len() * 24);
    for p in pts {
        out.extend_from_slice(&p.x.to_le_bytes());
        out.extend_from_slice(&p.y.to_le_bytes());
        out.extend_from_slice(&p.z.to_le_bytes());
    }
    out
}

fn sort_key(p: &Point3<f64>) -> (f64, f64, f64) {
    (p.x, p.y, p.z)
}
