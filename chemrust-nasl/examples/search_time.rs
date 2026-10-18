//! Runtime benchmark, part 2 of `docs/benchmarks.md`.
//!
//! Usage:
//!   cargo run -p chemrust-nasl --release --example search_time -- \
//!       <model.cell|model.cif> <bondlength> [--runs N] [--out FILE.json]

mod bench_common;

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
    for _ in 0..runs {
        let t = Instant::now();
        config.search_sites();
        runs_ms.push(t.elapsed().as_secs_f64() * 1e3);
    }

    let json = ResultJson {
        benchmark: "search_time",
        model: model_path.to_string_lossy().to_string(),
        bondlength: r,
        atom_count: points.len(),
        runs_ms,
        site_counts: counts,
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
