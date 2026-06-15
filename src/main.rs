//! Binary entry point: parse command-line arguments and hand work to the library.
use clap::Parser;
use kamino_cli::{run_with_args, Args};

fn main() -> anyhow::Result<()> {
    // Keep the binary thin so tests can exercise the same logic through the library.
    let args = Args::parse();
    run_with_args(args)
}
