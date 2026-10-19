//! Grid-search completeness baseline and accuracy-cost sweep,
//! part 1 of `docs/benchmarks.md`.
//!
//! A voxel grid fills the Cartesian box of the model atoms. A grid
//! point is a candidate when no atom is closer than r and at least
//! one atom sits in the window [r, r + delta]. The per-step metrics
//! compare the grid against the exact `SearchConfig` result.
//!
//! `--sweep N` runs the baseline at step, step/2, ..., step/2^(N-1).
//! The sweep exposes the trade-off: a coarse step misses exact sites
//! (incomplete), a fine step covers them but pays in candidate count
//! and scan time. The exact engine covers all sites in one
//! closed-form pass with no step parameter.
//!
//! Usage:
//!   cargo run -p chemrust-nasl --release --example grid_baseline -- \
//!       <model.cell|model.cif> <bondlength> [--step 0.05] [--delta 0.1] \
//!       [--sweep N] [--out FILE.json]

mod bench_common;

use std::fs::write;
use std::path::PathBuf;
use std::time::Instant;

use chemrust_nasl::{SearchConfig, SearchReports};
use kd_tree::KdIndexTree;
use nalgebra::{distance_squared, Point3, Vector3};
use rayon::prelude::*;
use serde::Serialize;

/// Tolerance on the squared-distance scale, matching the engine.
const TOL: f64 = 1.0e-5;

/// One candidate cell: grid indices plus the bond class.
/// Class 1 = single, 2 = double, 3 = multi.
type Cell = (usize, usize, usize, u8);

#[derive(Serialize)]
struct Counts {
    multi: usize,
    single: usize,
    double: usize,
}

/// Metrics for one grid step.
#[derive(Serialize)]
struct StepBlock {
    step: f64,
    scan_ms: f64,
    candidate_total: usize,
    candidate_counts: Counts,
    matched_sites: usize,
    missed_sites: usize,
    total_sites: usize,
    match_rate: f64,
    max_site_candidate_err: f64,
}

#[derive(Serialize)]
struct ResultJson {
    benchmark: &'static str,
    model: String,
    bondlength: f64,
    delta: f64,
    atom_count: usize,
    site_counts: Counts,
    steps: Vec<StepBlock>,
    context: bench_common::hardware_spec::Context,
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mut model = None;
    let mut bondlength: Option<f64> = None;
    let mut step = 0.05_f64;
    let mut delta = 0.1_f64;
    let mut sweep = 1;
    let mut out: Option<PathBuf> = None;
    let mut it = args.into_iter().skip(1);
    while let Some(a) = it.next() {
        match a.as_str() {
            "--step" => step = it.next().expect("--step takes a value").parse().expect("number"),
            "--delta" => delta = it.next().expect("--delta takes a value").parse().expect("number"),
            "--sweep" => sweep = it.next().expect("--sweep takes a value").parse().expect("number"),
            "--out" => out = Some(PathBuf::from(it.next().expect("--out takes a value"))),
            _ if model.is_none() => model = Some(PathBuf::from(a)),
            _ if bondlength.is_none() => {
                bondlength = Some(a.parse().expect("bondlength in Angstrom"))
            }
            _ => eprintln!("unexpected argument {a:?}"),
        }
    }
    let model_path = model.expect("usage: <model> <bondlength> [--step S] [--delta D] [--sweep N]");
    let r = bondlength.expect("usage: <model> <bondlength> [--step S] [--delta D] [--sweep N]");

    let structure = bench_common::bench_model::load_structure(&model_path).expect("load model");
    let points = bench_common::bench_model::cart_points(&structure);
    let to_check: Vec<(usize, Point3<f64>)> =
        points.iter().enumerate().map(|(i, p)| (i, *p)).collect();
    let config = SearchConfig::new(&to_check, &points, r);
    let report = config.search_sites();
    let sites = site_points(&report);

    let steps: Vec<f64> = (0..sweep.max(1)).map(|i| step / 2_f64.powi(i)).collect();
    let mut blocks = Vec::with_capacity(steps.len());
    for s in &steps {
        blocks.push(run_step(&points, &sites, r, delta, *s));
    }

    let json = ResultJson {
        benchmark: "grid_baseline",
        model: model_path.to_string_lossy().to_string(),
        bondlength: r,
        delta,
        atom_count: points.len(),
        site_counts: class_counts_from_report(&report),
        steps: blocks,
        context: bench_common::hardware_spec::Context::capture(rayon::current_num_threads()),
    };
    let text = serde_json::to_string_pretty(&json).expect("serialize result");
    match out {
        Some(path) => write(path, text).expect("write result file"),
        None => println!("{text}"),
    }
}

/// One grid step: scan and match against the exact sites.
fn run_step(points: &[Point3<f64>], sites: &[Point3<f64>], r: f64, delta: f64, step: f64) -> StepBlock {
    let t0 = Instant::now();
    let (cand_cells, cand_points) = scan_candidates(points, r, delta, step);
    let scan_ms = t0.elapsed().as_secs_f64() * 1e3;
    let cand_counts = class_counts(&cand_cells);
    let (matched, max_err) = match_sites(sites, &cand_points, delta);
    // `step` is intentionally not the match tolerance: completeness is
    // measured against a fixed physical tolerance, independent of the grid step.
    let total = sites.len();
    StepBlock {
        step,
        scan_ms,
        candidate_total: cand_cells.len(),
        candidate_counts: cand_counts,
        matched_sites: matched,
        missed_sites: total - matched,
        total_sites: total,
        match_rate: if total == 0 { 0.0 } else { matched as f64 / total as f64 },
        max_site_candidate_err: max_err,
    }
}

/// All exact sites of the report as plain points.
fn site_points(report: &SearchReports) -> Vec<Point3<f64>> {
    let mut sites = report
        .points()
        .cloned()
        .unwrap_or_default()
        .into_iter()
        .map(|p| p.point())
        .collect::<Vec<Point3<f64>>>();
    if let Some(v) = report.viable_single_points() {
        sites.extend(v.iter().map(|p| p.point()));
    }
    if let Some(v) = report.viable_double_points() {
        sites.extend(v.iter().map(|p| p.point()));
    }
    sites
}

/// Scan the voxel grid and return candidate cells plus their points.
fn scan_candidates(points: &[Point3<f64>], r: f64, delta: f64, step: f64) -> (Vec<Cell>, Vec<Point3<f64>>) {
    let mut lo = points[0];
    let mut hi = points[0];
    for p in points {
        lo.x = lo.x.min(p.x);
        lo.y = lo.y.min(p.y);
        lo.z = lo.z.min(p.z);
        hi.x = hi.x.max(p.x);
        hi.y = hi.y.max(p.y);
        hi.z = hi.z.max(p.z);
    }
    let grow = r + delta;
    lo -= Vector3::new(grow, grow, grow);
    hi += Vector3::new(grow, grow, grow);
    let dims = [
        ((hi.x - lo.x) / step).ceil().max(1.0) as usize,
        ((hi.y - lo.y) / step).ceil().max(1.0) as usize,
        ((hi.z - lo.z) / step).ceil().max(1.0) as usize,
    ];

    let atom_tree = KdIndexTree::build_by_ordered_float(points);
    let r2 = r * r;
    let win2 = (r + delta) * (r + delta);

    let rows: Vec<Vec<Cell>> = (0..dims[0])
        .into_par_iter()
        .map(|i| {
            let mut local: Vec<Cell> = Vec::new();
            for j in 0..dims[1] {
                for k in 0..dims[2] {
                    let p = Point3::new(
                        lo.x + i as f64 * step,
                        lo.y + j as f64 * step,
                        lo.z + k as f64 * step,
                    );
                    let mut in_window = 0;
                    let mut closer = 0;
                    for &id in atom_tree.within_radius(&p, r + delta) {
                        let d2 = distance_squared(&p, atom_tree.item(id));
                        if d2 < r2 - TOL {
                            closer += 1;
                        } else if d2 <= win2 + TOL {
                            in_window += 1;
                        }
                    }
                    if closer == 0 && in_window > 0 {
                        let class = if in_window == 1 {
                            1
                        } else if in_window == 2 {
                            2
                        } else {
                            3
                        };
                        local.push((i, j, k, class));
                    }
                }
            }
            local
        })
        .collect();
    let cand_cells: Vec<Cell> = rows.into_iter().flatten().collect();
    let cand_points = cand_cells
        .iter()
        .map(|&(i, j, k, _)| {
            Point3::new(
                lo.x + i as f64 * step,
                lo.y + j as f64 * step,
                lo.z + k as f64 * step,
            )
        })
        .collect();
    (cand_cells, cand_points)
}

/// Match exact sites to candidates. A site is matched when a
/// candidate sits within the fixed physical tolerance `tol`. That
/// tolerance is independent of the grid step. Returns the matched
/// count and the largest site-to-candidate error over matched sites.
fn match_sites(
    sites: &[Point3<f64>],
    cand_points: &[Point3<f64>],
    tol: f64,
) -> (usize, f64) {
    if sites.is_empty() {
        return (0, 0.0);
    }
    let cand_tree = KdIndexTree::build_by_ordered_float(cand_points);
    let mut matched = 0;
    let mut max_err = 0.0_f64;
    for p in sites {
        let hits = cand_tree.within_radius(p, tol);
        if hits.is_empty() {
            continue;
        }
        matched += 1;
        let err = hits
            .iter()
            .map(|&&id| distance_squared(&cand_points[id], p))
            .fold(f64::INFINITY, f64::min)
            .sqrt();
        max_err = max_err.max(err);
    }
    (matched, max_err)
}

fn class_counts(cells: &[Cell]) -> Counts {
    let mut c = Counts {
        multi: 0,
        single: 0,
        double: 0,
    };
    for cell in cells {
        match cell.3 {
            1 => c.single += 1,
            2 => c.double += 1,
            _ => c.multi += 1,
        }
    }
    c
}

fn class_counts_from_report(report: &SearchReports) -> Counts {
    Counts {
        multi: report.points().map(|v| v.len()).unwrap_or(0),
        single: report.viable_single_points().map(|v| v.len()).unwrap_or(0),
        double: report.viable_double_points().map(|v| v.len()).unwrap_or(0),
    }
}
