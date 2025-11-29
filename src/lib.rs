use clap::{ArgAction, Parser};

mod bubbles;
mod graph;
mod io;
mod output;
mod recode;
mod traverse;
//mod util;
mod variant_groups;

use recode::{DAYHOFF6_ALPHABET_SIZE, DAYHOFF6_BITS_PER_SYMBOL};

/// Build a node-based, colored de Bruijn graph from amino-acid proteomes and analyze bubbles.
#[derive(Parser, Debug)]
#[command(author, version, about, disable_version_flag = true)]
pub struct Args {
    /// Input directory with FASTA proteomes (plain or .gz)
    #[arg(short, long)]
    pub input: std::path::PathBuf,

    /// K-mer length on the Dayhoff6 alphabet [k=14]
    #[arg(short, long)]
    pub k: Option<usize>,

    /// Minimal sample frequency  [m=0.85]
    #[arg(short = 'm', long = "min-freq", default_value_t = 0.85, hide_default_value = true)]
    pub min_freq: f32,

    /// Maximum traversal depth from each start node [d=4]
    #[arg(short, long, default_value_t = 4, hide_default_value = true)]
    pub depth: usize,

    /// Output prefix [o=kamino]
    #[arg(short, long, default_value = "kamino")]
    pub output: std::path::PathBuf,

    /// Number of constant positions to keep from the in-bubble k-mer [c=3]
    #[arg(short, long)]
    pub constant: Option<usize>,

    /// Maximum number of middle positions per variant group [l=2*k]
    #[arg(short = 'l', long = "length-middle")]
    pub length_middle: Option<usize>,

    /// Number of threads [t=1]
    #[arg(short = 't', long)]
    pub threads: Option<usize>,

    /// Display version information.
    #[arg(short = 'v', long = "version", action = ArgAction::Version)]
    pub version: (),
}

pub fn run_with_args(args: Args) -> anyhow::Result<()> {
    // Fixed Dayhoff6 recoding scheme
    let sym_bits = DAYHOFF6_BITS_PER_SYMBOL;
    let alphabet_size = DAYHOFF6_ALPHABET_SIZE;

    // k defaults and limits for Dayhoff6
    let default_k = 14usize;
    let max_k = 21usize;
    let default_constant = 3usize;

    anyhow::ensure!(
        args.min_freq >= 0.6,
        "min_freq ({}) must be ≥ 0.6.",
        args.min_freq
    );
    let min_freq = args.min_freq;

    let k = args.k.unwrap_or(default_k);
    anyhow::ensure!(
        (2..=max_k).contains(&k),
        "k={} is invalid for {}: allowed range is 2..={} (default {})",
        k,
        "Dayhoff6",
        max_k,
        default_k
    );

    // Validate constant against k-1
    let k1 = k - 1;
    let constant = args.constant.unwrap_or_else(|| default_constant.min(k1));
    anyhow::ensure!(
        constant <= k1,
        "constant ({}) must be ≤ k-1 ({}).",
        constant,
        k1
    );

    let length_middle = args.length_middle.unwrap_or_else(|| 2 * k);
    anyhow::ensure!(
        length_middle >= 1,
        "length_middle ({}) must be ≥ 1.",
        length_middle
    );

    // Decide how many threads to use (default: 1)
    let num_threads: usize = args.threads.unwrap_or(1);
    anyhow::ensure!(num_threads >= 1, "threads must be ≥ 1");

    eprintln!("kamino v{}", env!("CARGO_PKG_VERSION"));
    eprintln!(
        "parameters: k={} constant={} min_freq={} depth={} length_middle={} threads={} input={} output={}",
        k,
        constant,
        min_freq,
        args.depth,
        length_middle,
        num_threads,
        args.input.display(),
        args.output.display()
    );

    // Build global graph
    let mut g = graph::Graph::new(k, sym_bits, alphabet_size);
    io::build_graph_from_dir(&args.input, k, &mut g, num_threads)?;

    // Basic stats
    io::print_graph_size(&g);
    //io::print_outdegree_histogram(&g);

    // Bubble endpoints
    let (start_kmers, end_kmers) = bubbles::find_bubble_endpoints(&g, min_freq, num_threads);
    eprintln!(
        "bubble endpoints: start={} end={}",
        start_kmers.len(),
        end_kmers.len()
    );

    // Traverse VGs
    let groups = traverse::find_variant_groups(
        &g,
        &start_kmers,
        &end_kmers,
        args.depth,
        min_freq,
        num_threads,
    );
    let total_paths: usize = groups.values().map(|v| v.len()).sum();
    eprintln!(
        "variant groups: groups={} paths={}",
        groups.len(),
        total_paths
    );

    // Emit amino-acid sequences per path for each variant group (call disabled, keep for debugging)
    // traverse::print_variant_group_sequences(&g, &groups);

    let (alignment_len, alignment_missing_pct) = output::write_outputs_with_head(
        &args.input,
        &args.output,
        &g,
        &groups,
        k,
        constant,
        min_freq,
        length_middle,
        num_threads,
    )?;
    let (fas_path, tsv_path, partitions_path) = output::output_paths(&args.output);
    eprintln!(
        "alignment: length={} missing={:.1}%",
        alignment_len, alignment_missing_pct
    );
    eprintln!(
        "output files:  {}, {}, and {}",
        fas_path.display(),
        tsv_path.display(),
        partitions_path.display()
    );

    Ok(())
}
