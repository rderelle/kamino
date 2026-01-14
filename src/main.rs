use clap::Parser;
use kamino_cli::{run_with_args, Args};

fn main() -> anyhow::Result<()> {
    let args = Args::parse();
    run_with_args(args)
}
