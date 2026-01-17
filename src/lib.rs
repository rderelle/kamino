//! # Kamino
//!
//! Kamino builds an amino-acid alignment in a reference-free, alignment-free manner from
//! a set of proteomes. It is not “better” than traditional marker-based pipelines, but it is
//! simpler and faster to use.
//!
//! Typical usage ranges from between-species to within-phylum phylogenetic analyses (bacteria,
//! archaea, and eukaryotes).
//!
//! ## Input modes
//! Kamino accepts proteome files as input in one of two modes:
//! - **Directory mode** (`--input-directory`): a directory containing FASTA proteomes
//!   (plain text or `.gz` compressed). Each file represents one isolate. Filenames minus the
//!   extension become sequence names in the final amino-acid alignment.
//! - **Table mode** (`--input-file`): a tab-delimited file mapping a species/sample name
//!   to a proteome path (one name + path pair per line). This is useful when file names
//!   do not encode the sample name or when proteomes are located in multiple directories.
//!
//!
//! ## Arguments
//! - `-k`, `--k`: k-mer length (default: 14; must be within the valid range for the
//!   selected recoding scheme).
//! - `-f`, `--min-freq`: minimum fraction of samples with an amino-acid per position
//!   (default: 0.85; must be ≥ 0.6).
//! - `-d`, `--depth`: maximum traversal depth from each start node (default: 6).
//! - `-o`, `--output`: output prefix for generated files (default: `kamino`).
//! - `-c`, `--constant`: number of constant positions retained from in-bubble k-mers
//!   (default: 3; must be ≤ k-1).
//! - `-l`, `--length-middle`: maximum number of middle positions per variant group
//!   (default: 2 * k; must be ≥ 1).
//! - `-m`, `--mask`: mask middle segments with long mismatch runs (default: 5).
//! - `-t`, `--threads`: number of threads used for graph construction and analysis
//!   (default: 1).
//! - `-r`, `--recode`: amino-acid recoding scheme (default: `sr6`).
//! - `-v`, `--version`: print version information and exit.
//!
//!
//! ## Important things to optimize
//! The main parameters governing the number of phylogenetic positions in the final alignment are
//! the k-mer size (-k), the depth of the recursive graph traversal (-d), and the minimum sample
//! frequency (-f). If the final alignment is too short, you may want to modify some
//! of these three parameters.
//!
//! Increasing the k-mer size may yield longer alignments, but only up to a point: bubble start and 
//! end k-mers shared by most samples become less frequent, resulting in fewer variant groups and a 
//! shorter final alignment.
//!
//! Increasing the depth of the recursive graph traversal (e.g., from 6 to 8) can also increase the size
//! of the final alignment as kamino detects more variant groups during the graph traversal. 
//!
//! Finally, it is also possible to produce larger alignments by decreasing the minimum fraction of samples
//! with an amino-acid (e.g., from 0.85 to 0.8), although samples will have more missing data in the final alignment.
//! Missing data in the alignment are represented by '-' (missing amino acid) and 'X' (ambiguous or masked amino acid).
//!
//!
//! ## Less important parameters
//! Besides testing/benchmarking, I would not recommend modifying these parameter values.
//!
//! The number of constant positions in the final alignment can be adjusted with the --constant parameter. These are
//! taken from the left flank of the end k-mer in each variant group, next to the middle positions. Because these
//! positions are recoded, some may become polymorphic once converted back to amino acids. Using the default c=3,
//! constant positions represent 60 to 75% of the alignment.
//!
//! The --mask parameter controls the amino-acid masking performed by kamino to prevent long runs of polymorphism from being
//! retained in the final alignment. These correspond to genuine but unwanted polymorphisms (e.g., micro-inversions) or,
//! less frequently, errors made by kamino (“misaligned” paths due to the presence of two consecutive indels). The minimum length
//! of polymorphism runs to be masked can be decreased using this parameter to be more stringent.
//!
//! The --length-middle parameter is used to filter out long variant groups. Increase this parameter to allow more
//! variant groups to be retained in the final alignment. Please note that this parameter might be removed in future versions.
//!
//! Finally, the 6-letter recoding scheme can be modified using the --recode parameter, although the default sr6 recoding
//! scheme performed the best in most of my tests (sr6 ≥ dayhoff6 ≫ kgb6).
//!
//!
//! ## Output files
//! The names of the three output files are controlled by a prefix (-o; default: `kamino`). The prefix
//! may include a directory path (e.g. `-o my_analyses/taxon1`). Note that the output directory is not
//! created by kamino and must already exist.
//!
//! The three output files are:
//! - `<prefix>_alignment.fas`: FASTA amino acid alignment of all samples.
//! - `<prefix>_missing.tsv`: Tab-delimited per-sample missingness percentages.
//! - `<prefix>_partitions.tsv`: Tab-delimited variant group coordinates (0-based) in the FASTA
//!   alignment, along with consensus protein names when the input proteomes are annotated.
//!
use clap::{ArgAction, Parser};

mod bubbles;
mod filter_groups;
mod graph;
mod io;
mod middle_mask;
mod output;
mod recode;
mod revert_aminoacid;
mod traverse;
//mod util;

pub use recode::RecodeScheme;

/// Build a node-based, colored de Bruijn graph from amino-acid proteomes and analyze bubbles.
#[derive(Parser, Debug)]
#[command(author, version, about, disable_version_flag = true)]
#[command(
    group = clap::ArgGroup::new("input_source")
        .required(true)
        .args(["input", "input_file"])
)]
pub struct Args {
    /// Input directory with FASTA proteomes (plain or .gz)
    #[arg(short, long = "input-directory")]
    pub input: Option<std::path::PathBuf>,

    /// Tab-delimited file mapping species name to proteome path
    #[arg(short = 'I', long = "input-file")]
    pub input_file: Option<std::path::PathBuf>,

    /// K-mer length [k=14]
    #[arg(short, long)]
    pub k: Option<usize>,

    /// Minimal fraction of samples with an amino-acid per position [m=0.85]
    #[arg(
        short = 'f',
        long = "min-freq",
        default_value_t = 0.85,
        hide_default_value = true
    )]
    pub min_freq: f32,

    /// Maximum traversal depth from each start node [d=6]
    #[arg(short, long, default_value_t = 6, hide_default_value = true)]
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

    /// Mask middle segments with long mismatch runs [m=5]
    #[arg(short = 'm', long = "mask", default_value_t = 5, hide_default_value = true)]
    pub mask: usize,

    /// Number of threads [t=1]
    #[arg(short = 't', long)]
    pub threads: Option<usize>,

    /// Recoding scheme [r=sr6]
    #[arg(short = 'r', long = "recode", value_enum, default_value_t = RecodeScheme::SR6, hide_default_value = true)]
    pub recode: RecodeScheme,

    /// Display version information.
    #[arg(short = 'v', long = "version", action = ArgAction::Version)]
    pub version: (),
}

pub fn run_with_args(args: Args) -> anyhow::Result<()> {
    // k defaults and limits for SR6
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
        args.recode,
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
    let (input_label, species_inputs) = if let Some(table) = args.input_file.as_ref() {
        (
            format!("input_file={}", table.display()),
            io::collect_species_inputs_from_table(table)?,
        )
    } else if let Some(dir) = args.input.as_ref() {
        (
            format!("input={}", dir.display()),
            io::collect_species_inputs_from_dir(dir)?,
        )
    } else {
        anyhow::bail!("Either --input-directory or --input-file must be provided.");
    };

    eprintln!(
        "parameters: k={} constant={} min_freq={} depth={} length_middle={} mask={} threads={} recode={} {} output={}",
        k,
        constant,
        min_freq,
        args.depth,
        length_middle,
        args.mask,
        num_threads,
        args.recode,
        input_label,
        args.output.display()
    );

    // Build global graph
    let mut g = graph::Graph::new(k, args.recode);
    io::build_graph_from_inputs(&species_inputs, k, &mut g, num_threads)?;

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
        &species_inputs,
        &args.output,
        &g,
        &groups,
        k,
        constant,
        min_freq,
        length_middle,
        args.mask,
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
