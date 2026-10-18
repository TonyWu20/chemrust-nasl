#![allow(dead_code)]
//! Host and build context written into every benchmark result file.
//!
//! See `docs/benchmarks.md`, part 2, for the field list.

use serde::Serialize;

/// One context block. Capture it next to the timed region.
#[derive(Debug, Clone, Serialize)]
pub struct Context {
    pub cpu_model: String,
    pub cpu_threads: usize,
    pub rayon_threads: usize,
    pub ram_mib: usize,
    pub os: String,
    pub kernel: String,
    pub rustc_version: String,
    pub chemrust_nasl_version: String,
    pub git_commit: String,
    pub freq_governor: String,
}

impl Context {
    /// Build the block. Pass the active rayon pool size.
    pub fn capture(rayon_threads: usize) -> Self {
        Self {
            cpu_model: model_name(),
            cpu_threads: std::thread::available_parallelism()
                .map(|n| n.get())
                .unwrap_or(0),
            rayon_threads,
            ram_mib: mem_total_mib(),
            os: std::env::consts::OS.to_string(),
            kernel: uname_kernel(),
            rustc_version: run_cmd("rustc", &["--version"]),
            chemrust_nasl_version: env!("CARGO_PKG_VERSION").to_string(),
            git_commit: run_cmd("git", &["rev-parse", "--short", "HEAD"]),
            freq_governor: read_file_trim("/sys/devices/system/cpu/cpu0/cpufreq/scaling_governor"),
        }
    }
}

fn model_name() -> String {
    if let Ok(text) = std::fs::read_to_string("/proc/cpuinfo") {
        for line in text.lines() {
            if let Some(rest) = line.strip_prefix("model name") {
                return rest.trim_start().trim_start_matches(':').trim().to_string();
            }
        }
    }
    "unknown".to_string()
}

fn mem_total_mib() -> usize {
    if let Ok(text) = std::fs::read_to_string("/proc/meminfo") {
        for line in text.lines() {
            if let Some(rest) = line.strip_prefix("MemTotal:") {
                let kb: f64 = rest
                    .trim()
                    .split_whitespace()
                    .next()
                    .unwrap_or("0")
                    .parse()
                    .unwrap_or(0.0);
                return (kb / 1024.0) as usize;
            }
        }
    }
    0
}

fn uname_kernel() -> String {
    if std::env::consts::OS == "linux" {
        return run_cmd("uname", &["-r"]);
    }
    "unknown".to_string()
}

fn run_cmd(bin: &str, args: &[&str]) -> String {
    match std::process::Command::new(bin).args(args).output() {
        Ok(out) if out.status.success() => String::from_utf8_lossy(&out.stdout).trim().to_string(),
        _ => "n/a".to_string(),
    }
}

fn read_file_trim(path: &str) -> String {
    std::fs::read_to_string(path)
        .map(|t| t.trim().to_string())
        .unwrap_or_else(|_| "not detected".to_string())
}
