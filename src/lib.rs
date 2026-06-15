//! # kamino
//!
//! kamino builds an amino-acid alignment in a reference-free, alignment-free manner from
//! a set of proteomes. It is not “better” than traditional marker-based pipelines, but it is
//! simpler and faster to use.
//!
//! Typical usage ranges from between-species to within-phylum phylogenetic analyses (bacteria,
//! archaea, and eukaryotes).
//!
//! ## Input modes
//! kamino accepts proteome files as input in one of two modes:
//! - **Directory mode** (`--input-directory`): a directory containing FASTA proteomes
//!   (plain text or `.gz` compressed). Each file represents one isolate. Filenames minus the
//!   extension become sequence names in the final amino-acid alignment.
//! - **Table mode** (`--input-file`): a tab-delimited file mapping a species/sample name
//!   to a proteome path (one name + path pair per line). This is useful when file names
//!   do not encode the sample name or when proteomes are located in multiple directories.
//!
//! In the directory mode, files are recognized by their extension (.fas, .fasta, .faa, .fa, .fna; gzipped ot not).
//!
//! For **bacterial** isolates, the phylogenomic alignment can also be generated directly from genome assemblies
//! by selecting the option `--genomes` (using either `-i` or `-I`). In this case, an ultra-fast but approximate
//! protein prediction is performed, and the predicted proteomes are written to a temporary directory.
//!
//! ## Arguments
//! - `-i`, --input-directory: input directory with FASTA proteomes (plain or .gz)
//! - `-I`, --input-file: tab-delimited file mapping species name to proteome path
//! - `--genomes`: treat input files as bacterial genomes and predict proteomes before analysis
//! - `-k`, `--k`: k-mer length [k=13]
//! - `-f`, `--min-freq`: minimum fraction of samples with an amino-acid per position [f=0.85]
//! - `-o`, `--output`: output prefix for generated files [o=`kamino`]
//! - `-c`, `--constant`: number of constant positions retained from in-bubble k-mers [c=3]
//! - `-l`, `--length-middle`: maximum number of middle positions per variant group [l=35]
//! - `-m`, `--mask`: mask middle segments with long mismatch runs [m=5]
//! - `-t`, `--threads`: number of threads [t=1]
//! - `-r`, `--recode`: amino-acid recoding scheme [r=`sr6`]
//! - `--nj`: generate a NJ tree from kamino alignment [nj=false]
//! - `-v`, `--version`: print version information and exit.
//!
//!
//! ## Optimising alignment size
//!
//! The final alignment size can mainly be increased in two ways:
//!
//! 1. Decrease the minimum fraction of samples required to carry an amino acid with
//!    `--min-freq`, for example from 0.85 to 0.80. This will produce larger
//!    alignments, but at the cost of increased missing data. Missing data are
//!    represented by `-` for missing amino acids and `X` for ambiguous or masked
//!    amino acids.
//!
//! 2. Increase the `--length-middle` parameter, which controls the maximum number of
//!    middle positions in variant groups, for example from 35 to 70. This allows longer
//!    variant groups to be retained in the final alignment.
//!
//! Conversely, if the alignment is too large, you can increase the minimum fraction
//! of samples and/or reduce the maximum length of middle positions.
//!
//!
//! ## Less important parameters
//!
//! Except for testing and benchmarking, I do not recommend changing these parameter
//! values.
//!
//! The default k-mer size has been chosen to maximise the final alignment length in most
//! conditions. Increasing it usually does not substantially increase the number of variant
//! groups.
//!
//! The number of constant positions in the final alignment can be adjusted with the
//! --constant parameter. These positions are taken from the left flank of the end
//! k-mer in each variant group, next to the middle positions. Because these positions
//! are recoded, some may become polymorphic once converted back to amino acids. With
//! the default value of c = 3, constant positions represent about 50% of the
//! alignment.
//!
//! The --mask parameter controls the amino-acid masking performed by kamino to
//! prevent long runs of polymorphism from being retained in the final alignment.
//! These runs correspond to genuine but unwanted polymorphisms, such as
//! micro-inversions, or, less frequently, errors made by kamino, such as misaligned
//! paths caused by two consecutive indels. The minimum length of polymorphic runs to
//! be masked can be decreased with this parameter to make the filtering more
//! stringent.
//!
//! Finally, the 6-letter recoding scheme can be modified with the --recode
//! parameter, although the default sr6 recoding scheme performed best in most of my
//! tests (sr6 >= dayhoff6 >> kgb6).
//!
//!
//! ## Output files
//! The names of the output files are controlled by a prefix (-o; default: `kamino`). The prefix
//! may include a directory path (e.g. `-o my_analyses/taxon1`). Note that the output directory is not
//! created by kamino and must already exist.
//!
//! The three output files are:
//! - `<prefix>_alignment.fas`: FASTA amino acid alignment of all samples.
//! - `<prefix>_missing.tsv`: Tab-delimited per-sample missingness percentages.
//! - `<prefix>_partitions.tsv`: Tab-delimited variant group coordinates (0-based) in the FASTA
//!   alignment, along with consensus protein names when the input proteomes are annotated.
//!
//! Additionally, a Neighbor-Joining (NJ) tree can be produced from the amino acid alignment
//! when the `--nj` argument is specified. Pairwise distances are computed using an F81
//! correction with LG stationary amino-acid frequencies. The resulting tree provides an
//! overview of isolate relationships and is not intended for detailed phylogenetic inference.
//!
use anyhow::Context;
use clap::Parser;

mod group_extraction;
mod group_filtering;
mod group_sorting;
mod io;
mod output;
mod phylo;
mod proba_filter;
mod protein_prediction;
mod recode;

pub use recode::RecodeScheme;

#[derive(Parser, Debug)]
#[command(name = "kamino", author, version, about)]
#[command(group = clap::ArgGroup::new("input_source").required(true).multiple(true).args(["input", "input_file"]))]
/// Parsed command-line options shared by the binary and integration tests.
pub struct Args {
    /// Directory containing one proteome FASTA per species.
    #[arg(short, long = "input-directory")]
    pub input: Option<std::path::PathBuf>,
    /// TSV table with `species_name<TAB>path_to_fasta` rows.
    #[arg(short = 'I', long = "input-file")]
    pub input_file: Option<std::path::PathBuf>,
    /// Treat inputs as nucleotide genomes and predict proteins first.
    #[arg(long = "genomes")]
    pub genomes: bool,
    /// k-mer size used for anchor extraction [k=13].
    #[arg(short, long)]
    pub k: Option<usize>,
    /// Minimum fraction of species required to support an anchor or output column [f=0.85].
    #[arg(
        short = 'f',
        long = "min-freq",
        default_value_t = 0.85,
        hide_default_value = true
    )]
    pub min_freq: f32,
    /// Output prefix used to derive alignment, missing-data, partition, and tree paths.
    #[arg(short, long, default_value = "kamino")]
    pub output: std::path::PathBuf,
    /// Number of constant right-anchor columns to append after each variable block [c=3].
    #[arg(short, long)]
    pub constant: Option<usize>,
    /// Maximum amino-acid length allowed between two adjacent shared anchors [l=35].
    #[arg(
        short = 'l',
        long = "length-middle",
        default_value_t = 35,
        hide_default_value = true
    )]
    pub length_middle: usize,
    /// Consecutive amino-acid differences from the group consensus required to mask a row [m=5]; 0 disables.
    #[arg(
        short = 'm',
        long = "mask",
        default_value_t = 5,
        hide_default_value = true
    )]
    pub mask: usize,
    /// Number of worker threads used by parallel stages [t=1].
    #[arg(short = 't', long, default_value_t = 1, hide_default_value = true)]
    pub threads: usize,
    /// Six-state amino-acid recoding scheme used before k-mer encoding [r=sr6].
    #[arg(short='r', long="recode", value_enum, default_value_t=RecodeScheme::SR6, hide_default_value=true)]
    pub recode: RecodeScheme,
    /// Also write a neighbor-joining tree inferred from the concatenated alignment.
    #[arg(long = "nj")]
    pub nj: bool,
}

fn print_startup_banner(args: &Args, k: usize, constant: usize) {
    let mut parameters = Vec::new();

    if let Some(input) = args.input.as_ref() {
        parameters.push(format!("input-directory={}", input.display()));
    }
    if let Some(input_file) = args.input_file.as_ref() {
        parameters.push(format!("input-file={}", input_file.display()));
    }
    if args.genomes {
        parameters.push("genomes=true".to_string());
    }
    parameters.push(format!("k={k}"));
    parameters.push(format!("min-freq={}", args.min_freq));
    parameters.push(format!("output={}", args.output.display()));
    parameters.push(format!("constant={constant}"));
    parameters.push(format!("length-middle={}", args.length_middle));
    parameters.push(format!("mask={}", args.mask));
    parameters.push(format!("threads={}", args.threads));
    parameters.push(format!("recode={}", args.recode));
    if args.nj {
        parameters.push("nj=true".to_string());
    }

    eprintln!("kamino {}", env!("CARGO_PKG_VERSION"));
    eprintln!("parameters: {}", parameters.join(" "));
}

/// Validate arguments, collect inputs, run the analysis, and write output files.
pub fn run_with_args(args: Args) -> anyhow::Result<()> {
    // Defaults and bounds are kept here so tests and the CLI share identical behavior.
    let default_k = 13usize;
    let max_k = 21usize;
    let k = args.k.unwrap_or(default_k);
    let constant = args.constant.unwrap_or(3usize.min(k));
    print_startup_banner(&args, k, constant);
    anyhow::ensure!(
        (0.6..=1.0).contains(&args.min_freq),
        "min_freq must be between 0.6 and 1.0"
    );
    anyhow::ensure!(args.threads > 0, "threads must be >=1");
    anyhow::ensure!((1..=max_k).contains(&k), "invalid k");
    anyhow::ensure!(constant <= k, "constant <= k");
    // Merge the optional input sources into one sorted list of species inputs.
    let mut species_inputs = Vec::new();
    if let Some(t) = args.input_file.as_ref() {
        species_inputs.extend(io::collect_species_inputs_from_table(t)?);
    }
    if let Some(d) = args.input.as_ref() {
        species_inputs.extend(io::collect_species_inputs_from_dir(d)?);
    }
    let output_dir = args
        .output
        .parent()
        .map(std::path::Path::to_path_buf)
        .unwrap_or_else(|| std::path::PathBuf::from("."));

    eprintln!("# process input files");
    // Genome mode creates temporary predicted proteomes next to the output prefix;
    // keeping the TempDir alive until the end guarantees those files remain readable.
    let mut genomes_tmpdir = None;
    if args.genomes {
        eprintln!(" . predict proteins from bacterial genomes");
        let tmpdir = tempfile::Builder::new()
            .prefix("tmp_kamino_")
            .rand_bytes(6)
            .tempdir_in(&output_dir)
            .with_context(|| format!("create temporary directory in {}", output_dir.display()))?;
        species_inputs =
            protein_prediction::predict_proteomes(&species_inputs, tmpdir.path(), args.threads)?;
        genomes_tmpdir = Some(tmpdir);
    }
    anyhow::ensure!(
        !species_inputs.is_empty(),
        "At least one input source with files must be provided."
    );
    // The extraction, sorting, and filtering stages return in-memory rows and partition metadata;
    // the output module owns all filesystem side effects after this point.
    let raw_groups = group_extraction::extract_groups(
        &species_inputs,
        k,
        args.min_freq,
        args.length_middle,
        args.recode,
        args.threads,
    )?;
    eprintln!("# analyse variant groups");
    eprintln!(" . raw variant groups: {}", raw_groups.groups.len());

    let sorted_groups = group_sorting::sort_and_deduplicate_groups(raw_groups, args.recode)?;
    eprintln!(
        " . sorted variant groups: {}",
        sorted_groups.raw_candidates.len()
    );

    let res = group_filtering::filter_groups(
        sorted_groups,
        args.min_freq,
        constant,
        args.mask,
        args.threads,
    )?;
    eprintln!(" . filtered variant groups: {}", res.partitions.len());

    let (alen, amiss) = output::write_outputs(
        &args.output,
        &res.species_names,
        res.concat,
        res.partitions,
        res.partition_names,
        args.nj,
        args.threads,
    )?;

    eprintln!("# output files");
    eprintln!(" . alignment: length={} missing={:.1}%", alen, amiss);
    if args.nj {
        eprintln!(" . NJ tree");
    }
    drop(genomes_tmpdir);
    Ok(())
}
