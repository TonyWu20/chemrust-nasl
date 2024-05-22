use clap::{Parser, ValueEnum};

#[derive(Parser)]
#[command(author, version, about, long_about= None)]
pub struct Args {
    #[arg(short, long)]
    pub mode: Option<ProgramMode>,
    pub config_loc: Option<String>,
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
/// Let the user decide to run with configs written in a prepared `yaml`,
/// or in interactive mode.
pub enum ProgramMode {
    /// Config
    C,
    /// Interactive
    I,
}
