//! Scalability on surface area: time the exact engine and a
//! z-pinned voxel-grid baseline across in-plane supercell sizes.
//!
//! The slab keeps a fixed layer count and vacuum gap. The z axis is
//! pinned to the surface + vacuum band, mirroring the app per-axis
//! range flags (`chemrust-nasl-app` `run_options.rs` / `search_job.rs`):
//! only the surface-band atoms act as bond anchors, and the grid
//! candidate space is restricted to the same z band. The scaling
//! variable is the in-plane supercell, i.e. the surface area.
//!
//! This isolates the quantity that matters for high-throughput
//! surface screening: cost should track surface area, not layer
//! depth.
//!
//! Usage:
//!   cargo run -p chemrust-nasl --release --example scalability -- \
//!       [--species Pt] [--a 3.924] <h> <k> <l> [--layers 4] \
//!       [--inplane 1,2,3,4] [--vacuum 12.0] [--bondlength 2.5] \
//!       [--delta 0.1] [--step 0.1] [--grid] [--surface] \
//!       [--zmin F] [--zmax F] [--runs 3] [--probe] [--out FILE.json]

mod bench_common;

use std::collections::{HashMap, VecDeque};
use std::fs::write;
use std::path::PathBuf;
use std::str::FromStr;
use std::time::Instant;

use castep_periodic_table::element::ElementSymbol;
use chemrust_geometry::{slab::fcc_bulk, transform::SurfaceRotation, Structure};
use chemrust_nasl::SearchConfig;
use kd_tree::KdIndexTree;
use nalgebra::{distance_squared, Point3};
use rayon::prelude::*;
use serde::Serialize;

/// Tolerance on the squared-distance scale, matching the engine.
const TOL: f64 = 1.0e-5;

/// One candidate cell: grid indices plus the bond class.
type Cell = (usize, usize, usize, u8);

#[derive(Serialize)]
struct GridBlock {
    step: f64,
    scan_ms: f64,
    candidates: usize,
    clusters: usize,
    matched_sites: usize,
    missed_sites: usize,
}

#[derive(Serialize)]
struct Row {
    in_plane: String,
    atom_count: usize,
    anchor_count: usize,
    runs_ms: Vec<f64>,
    multi: usize,
    single: usize,
    double: usize,
    grid: Option<GridBlock>,
}

#[derive(Serialize)]
struct ResultJson {
    benchmark: &'static str,
    species: String,
    hkl: [i32; 3],
    layers: usize,
    bondlength: f64,
    vacuum: f64,
    z_band: [f64; 2],
    rows: Vec<Row>,
    context: bench_common::hardware_spec::Context,
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mut species = "Pt".to_string();
    let mut a_latt = 3.924_f64;
    let mut hkl: [i32; 3] = [1, 1, 1];
    let mut hkl_it = 0;
    let mut layers = 4usize;
    let mut inplane = vec![1usize, 2, 3, 4];
    let mut vacuum = 12.0_f64;
    let mut bondlength = 2.5_f64;
    let mut delta = 0.1_f64;
    let mut step = 0.1_f64;
    let mut runs = 3usize;
    let mut grid = false;
    let mut surface = true;
    let mut zmin: Option<f64> = None;
    let mut zmax: Option<f64> = None;
    let mut out: Option<PathBuf> = None;
    let mut probe = false;
    let mut it = args.into_iter().skip(1);
    while let Some(a) = it.next() {
        match a.as_str() {
            "--species" => species = it.next().expect("--species takes a value"),
            "--a" => a_latt = it.next().expect("--a takes a value").parse().expect("number"),
            "--layers" => layers = it.next().expect("--layers takes a value").parse().expect("number"),
            "--inplane" => inplane = it
                .next()
                .expect("--inplane takes a value")
                .split(',')
                .map(|s| s.parse().expect("integer"))
                .collect(),
            "--vacuum" => vacuum = it.next().expect("--vacuum takes a value").parse().expect("number"),
            "--bondlength" => bondlength = it.next().expect("--bondlength takes a value").parse().expect("number"),
            "--delta" => delta = it.next().expect("--delta takes a value").parse().expect("number"),
            "--step" => step = it.next().expect("--step takes a value").parse().expect("number"),
            "--runs" => runs = it.next().expect("--runs takes a value").parse().expect("number"),
            "--zmin" => zmin = Some(it.next().expect("--zmin takes a value").parse().expect("number")),
            "--zmax" => zmax = Some(it.next().expect("--zmax takes a value").parse().expect("number")),
            "--grid" => grid = true,
            "--surface" => surface = true,
            "--no-surface" => surface = false,
            "--probe" => probe = true,
            "--out" => out = Some(PathBuf::from(it.next().expect("--out takes a value"))),
            _ if hkl_it < 3 => {
                hkl[hkl_it] = a.parse().expect("hkl index");
                hkl_it += 1;
            }
            _ => eprintln!("unexpected argument {a:?}"),
        }
    }
    // hkl is optional; [1,1,1] is the default when no positional is given.

    let sym = ElementSymbol::from_str(&species).expect("a valid element symbol for --species");
    if probe {
        let structure = slab(layers, 1, 1, &sym, a_latt, hkl, vacuum);
        let mut c: Vec<f64> = structure
            .frac_coords
            .iter()
            .map(|f| f.0.z)
            .collect();
        c.sort_by(|a, b| a.total_cmp(b));
        c.dedup_by(|a, b| (*a - *b).abs() < 1e-6);
        let c_len = structure.cell.as_ref().unwrap().tensor().column(2).norm();
        println!("c-lattice-len={:.4} A", c_len);
        for v in c {
            println!("plane c-frac={v:.4}");
        }
        std::process::exit(0);
    }

    let mut rows: Vec<Row> = Vec::with_capacity(inplane.len());
    let mut z_band_used: [f64; 2] = [0.0, 1.0];
    for s in &inplane {
        let structure = slab(layers, *s, *s, &sym, a_latt, hkl, vacuum);
        let points = bench_common::bench_model::cart_points(&structure);
        let frac = &structure.frac_coords;

        // Resolve the z band. An explicit --zmin/--zmax wins.
        // Otherwise, when --surface, pin to [top plane c-frac, 1.0].
        let top_frac = frac.iter().map(|f| f.0.z).fold(0.0_f64, f64::max);
        let (zmin, zmax) = if zmin.is_some() || zmax.is_some() {
            (zmin.unwrap_or(0.0), zmax.unwrap_or(1.0))
        } else if surface {
            (top_frac - 1.0e-3, 1.0)
        } else {
            (0.0, 1.0)
        };
        z_band_used = [zmin, zmax];

        // Anchors: atoms whose c-frac is inside the band.
        let to_check: Vec<(usize, Point3<f64>)> = points
            .iter()
            .enumerate()
            .filter(|(i, _)| {
                let cf = frac[*i].0.z;
                cf >= zmin - 1.0e-9 && cf <= zmax + 1.0e-9
            })
            .map(|(i, p)| (i, *p))
            .collect();

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
        let exact_sites = site_points(&report);
        let mean_ms = runs_ms.iter().sum::<f64>() / runs.max(1) as f64;

        let grid_block = if grid {
            let c_len = structure.cell.as_ref().unwrap().tensor().column(2).norm();
            let z_lo = zmin * c_len;
            let z_hi = zmax * c_len;
            let (cand_cells, cand_points, scan_ms) =
                z_pinned_scan(&points, bondlength, delta, step, z_lo, z_hi);
            let clusters = cluster_candidates(&cand_cells);
            let (matched, missed) = match_sites(&exact_sites, &cand_points, &cand_cells, delta);
            Some(GridBlock {
                step,
                scan_ms,
                candidates: cand_cells.len(),
                clusters,
                matched_sites: matched,
                missed_sites: missed,
            })
        } else {
            None
        };

        rows.push(Row {
            in_plane: format!("{s}x{s}"),
            atom_count: points.len(),
            anchor_count: to_check.len(),
            runs_ms,
            multi,
            single,
            double,
            grid: grid_block,
        });
        let g = rows.last().unwrap().grid.as_ref();
        match g {
            Some(g) => println!(
                "{species}{hkl:?} s={s} atoms={:3} anchors={:3}  \
                 multi/single/double={}/{}/{}  exact~{:.1} ms  \
                 grid(step={}) cands={} clusters={} matched={}/{} {}ms",
                points.len(),
                to_check.len(),
                multi,
                single,
                double,
                mean_ms,
                g.step,
                g.candidates,
                g.clusters,
                g.matched_sites,
                exact_sites.len(),
                g.scan_ms
            ),
            None => println!(
                "{species}{hkl:?} s={s} atoms={:3} anchors={:3}  \
                 multi/single/double={}/{}/{}  exact~{:.1} ms",
                points.len(),
                to_check.len(),
                multi,
                single,
                double,
                mean_ms
            ),
        }
    }

    let json = ResultJson {
        benchmark: "scalability",
        species: species.clone(),
        hkl,
        layers,
        bondlength,
        vacuum,
        z_band: z_band_used,
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

/// Voxel scan over the in-plane box with the z axis pinned to
/// `[z_lo, z_hi]` (Cartesian). A cell is a candidate when no atom
/// is closer than r and at least one atom sits in [r, r + delta].
/// Mirrors `grid_baseline::scan_candidates` with a pinned z range.
fn z_pinned_scan(
    points: &[Point3<f64>],
    r: f64,
    delta: f64,
    step: f64,
    z_lo: f64,
    z_hi: f64,
) -> (Vec<Cell>, Vec<Point3<f64>>, f64) {
    let grow = r + delta;
    let mut lo = points[0];
    let mut hi = points[0];
    for p in points {
        lo.x = lo.x.min(p.x);
        lo.y = lo.y.min(p.y);
    }
    for p in points {
        hi.x = hi.x.max(p.x);
        hi.y = hi.y.max(p.y);
    }
    lo.x -= grow;
    lo.y -= grow;
    lo.z = z_lo;
    hi.x += grow;
    hi.y += grow;
    hi.z = z_hi;
    let nx = ((hi.x - lo.x) / step).ceil().max(1.0) as usize;
    let ny = ((hi.y - lo.y) / step).ceil().max(1.0) as usize;
    let nz = ((hi.z - lo.z) / step).ceil().max(1.0) as usize;
    let tree = KdIndexTree::build_by_ordered_float(points);
    let r2 = r * r;
    let outer2 = (r + delta) * (r + delta);
    let t0 = Instant::now();
    let rows: Vec<Vec<Cell>> = (0..nx)
        .into_par_iter()
        .map(|i| {
            let mut local: Vec<Cell> = Vec::new();
            for j in 0..ny {
                for k in 0..nz {
                    let p = Point3::new(
                        lo.x + i as f64 * step,
                        lo.y + j as f64 * step,
                        lo.z + k as f64 * step,
                    );
                    let mut closer = 0;
                    let mut in_window = 0;
                    for &id in tree.within_radius(&p, r + delta) {
                        let d2 = distance_squared(&p, tree.item(id));
                        if d2 < r2 - TOL {
                            closer += 1;
                        } else if d2 <= outer2 + TOL {
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
    let scan_ms = t0.elapsed().as_secs_f64() * 1e3;
    let cand_points: Vec<Point3<f64>> = cand_cells
        .iter()
        .map(|&(i, j, k, _)| {
            Point3::new(
                lo.x + i as f64 * step,
                lo.y + j as f64 * step,
                lo.z + k as f64 * step,
            )
        })
        .collect();
    (cand_cells, cand_points, scan_ms)
}

/// Collect all exact site positions from a report.
fn site_points(report: &chemrust_nasl::SearchReports) -> Vec<Point3<f64>> {
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

/// 26-neighborhood cluster count over candidate cells.
fn cluster_candidates(cand_cells: &[Cell]) -> usize {
    if cand_cells.is_empty() {
        return 0;
    }
    let index: HashMap<(usize, usize, usize), usize> = cand_cells
        .iter()
        .enumerate()
        .map(|(n, &(i, j, k, _))| ((i, j, k), n))
        .collect();
    let mut visited = vec![false; cand_cells.len()];
    let mut clusters = 0;
    for start in 0..cand_cells.len() {
        if visited[start] {
            continue;
        }
        visited[start] = true;
        clusters += 1;
        let mut queue = VecDeque::new();
        queue.push_back(start);
        while let Some(cur) = queue.pop_front() {
            let (ci, cj, ck, _) = cand_cells[cur];
            for di in -1_i32..=1 {
                for dj in -1_i32..=1 {
                    for dk in -1_i32..=1 {
                        if di == 0 && dj == 0 && dk == 0 {
                            continue;
                        }
                        let ni = (ci as i32 + di).max(0) as usize;
                        let nj = (cj as i32 + dj).max(0) as usize;
                        let nk = (ck as i32 + dk).max(0) as usize;
                        if let Some(&n) = index.get(&(ni, nj, nk)) {
                            if !visited[n] {
                                visited[n] = true;
                                queue.push_back(n);
                            }
                        }
                    }
                }
            }
        }
    }
    clusters
}

/// Match exact sites to candidate clusters. A site is matched when a
/// candidate sits within the fixed physical tolerance `tol`. Returns
/// (matched, missed).
fn match_sites(
    sites: &[Point3<f64>],
    cand_points: &[Point3<f64>],
    _cand_cells: &[Cell],
    tol: f64,
) -> (usize, usize) {
    if cand_points.is_empty() {
        return (0, sites.len());
    }
    let cand_tree = KdIndexTree::build_by_ordered_float(cand_points);
    let mut matched = 0usize;
    for p in sites {
        if !cand_tree.within_radius(p, tol).is_empty() {
            matched += 1;
        }
    }
    (matched, sites.len() - matched)
}
