//! Build phylogenomic amino-acid alignments directly from proteomes.
//!
//! `kamino` finds exact k-mer anchors shared across samples, extracts the variable regions
//! between adjacent anchors, and concatenates the retained regions into an alignment. This
//! avoids choosing a reference genome or first defining a fixed set of marker genes. It is
//! intended as a fast way to generate datasets for comparisons ranging from closely related
//! isolates to broader bacterial, archaeal, and eukaryotic groups.
//!
//! This crate contains the command-line application's reusable entry point. Most users should
//! install and run the `kamino` executable; callers embedding it can construct [`Args`] and
//! pass them to [`run_with_args`].
//!
//! # Quick start
//!
//! Put one proteome FASTA file per sample in a directory, then run:
//!
//! ```text
//! kamino --input-directory proteomes --output results/run1 --threads 4
//! ```
//!
//! The output directory (`results` above) must already exist. Run `kamino --help` for the
//! complete and authoritative list of options and their defaults.
//!
//! # Input
//!
//! One or both of the following input sources may be supplied:
//!
//! - **Directory input** (`-i`, `--input-directory`) reads every `.fa`, `.fas`, `.fasta`,
//!   `.faa`, or `.fna` file in a directory. The same extensions followed by `.gz` are also
//!   accepted. A filename without its FASTA and compression extensions becomes the sample
//!   name (for example, `isolate.1.faa.gz` becomes `isolate.1`).
//! - **Manifest input** (`-I`, `--input-file`) reads a headerless, two-column TSV containing
//!   `sample_name<TAB>proteome_path` on each line. Relative paths are resolved from the
//!   manifest's directory, so a manifest can be moved together with its input files.
//!
//! For bacterial genome assemblies rather than predicted proteomes, add `--genomes`.
//! Kamino then performs a fast, approximate protein prediction in a temporary directory
//! before building the alignment.
//!
//! # Output
//!
//! `--output` is a filename prefix, not a directory. With the default prefix `kamino`, the
//! application writes:
//!
//! - `kamino_alignment.fas` — the concatenated amino-acid alignment;
//! - `kamino_missing.tsv` — the percentage of missing or masked sites per sample; and
//! - `kamino_partitions.tsv` — zero-based coordinates and consensus protein names for the
//!   retained variant groups.
//!
//! `--nj` additionally writes `kamino_NJ.tree`, a neighbor-joining tree in Newick format.
//! It uses F81-corrected distances with LG stationary amino-acid frequencies and is intended
//! as a quick overview, not as a replacement for detailed phylogenetic inference.
//!
//! # Choosing parameters
//!
//! The defaults are suitable for most analyses. The main control on alignment size is
//! `-f`/`--min-freq`: lower values retain sites present in fewer samples, producing a longer
//! alignment with more missing data. Missing residues are written as `-`; ambiguous or masked
//! residues are written as `X`.
//!
//! The other parameters refine how anchors and retained regions are handled:
//!
//! - `-k`/`--k` sets the anchor k-mer size. It can be adjusted for the evolutionary scale of
//!   the dataset; for example, `-k 9` can be useful for comparisons within a species.
//! - `-c`/`--constant` controls how many constant anchor sites are added to each retained
//!   variant group in the final alignment. It does not affect which variants are identified.
//! - `-l`/`--length-middle` sets the maximum amino-acid length allowed between adjacent shared
//!   anchors.
//! - `-m`/`--mask` sets the minimum run of consecutive amino-acid differences to mask; use `0`
//!   to disable this filtering.
//! - `-r`/`--recode` enables the Dayhoff6, SR6, or KGB6 six-state recoding scheme for anchor
//!   discovery. Output is still written using the original amino-acid alphabet.
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
    /// Directory containing proteome FASTA files.
    #[arg(short, long = "input-directory", help_heading = "Main parameters")]
    pub input: Option<std::path::PathBuf>,

    /// TSV table with `species_name<TAB>path_to_fasta` rows.
    #[arg(short = 'I', long = "input-file", help_heading = "Main parameters")]
    pub input_file: Option<std::path::PathBuf>,

    /// Prefix for output files
    #[arg(
        short,
        long,
        default_value = "kamino",
        help_heading = "Main parameters"
    )]
    pub output: std::path::PathBuf,

    /// k-mer size [default: 8; 13 when six-state recoding is enabled].
    #[arg(short, long, help_heading = "Main parameters")]
    pub k: Option<usize>,

    /// Optional six-state amino-acid k-mer recoding scheme.
    #[arg(
        short = 'r',
        long = "recode",
        value_enum,
        help_heading = "Main parameters"
    )]
    pub recode: Option<RecodeScheme>,

    /// Minimum fraction of species present per alignment position [f=0.85].
    #[arg(
        short = 'f',
        long = "min-freq",
        default_value_t = 0.85,
        hide_default_value = true,
        help_heading = "Main parameters"
    )]
    pub min_freq: f32,

    /// Constant positions [default: 1; 3 when six-state recoding is enabled].
    #[arg(short, long, help_heading = "Main parameters")]
    pub constant: Option<usize>,

    /// Maximum amino-acid length between two adjacent shared anchors [l=35].
    #[arg(
        short = 'l',
        long = "length-middle",
        default_value_t = 35,
        hide_default_value = true,
        help_heading = "Main parameters"
    )]
    pub length_middle: usize,

    /// Consecutive amino-acid differences required for masking [m=5]; 0 disables.
    #[arg(
        short = 'm',
        long = "mask",
        default_value_t = 5,
        hide_default_value = true,
        help_heading = "Main parameters"
    )]
    pub mask: usize,

    /// Number of threads [t=1].
    #[arg(
        short = 't',
        long,
        default_value_t = 1,
        hide_default_value = true,
        help_heading = "Main parameters"
    )]
    pub threads: usize,

    /// Treat inputs as bacterial genomes and predict proteins first.
    #[arg(long = "genomes", help_heading = "Optional input")]
    pub genomes: bool,

    /// Build a neighbor-joining tree.
    #[arg(long = "nj", help_heading = "Optional output")]
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
    parameters.push(format!(
        "recode={}",
        args.recode
            .map_or_else(|| "none".to_string(), |scheme| scheme.to_string())
    ));
    if args.nj {
        parameters.push("nj=true".to_string());
    }

    eprintln!("kamino {}", env!("CARGO_PKG_VERSION"));
    eprintln!("parameters: {}", parameters.join(" "));
}

/// Validate arguments, collect inputs, run the analysis, and write output files.
pub fn run_with_args(args: Args) -> anyhow::Result<()> {
    // Defaults and bounds are kept here so tests and the CLI share identical behavior.
    let (alphabet, k, constant) = resolve_alphabet_parameters(args.recode, args.k, args.constant)?;
    print_startup_banner(&args, k, constant);
    anyhow::ensure!(
        (0.6..=1.0).contains(&args.min_freq),
        "min_freq must be between 0.6 and 1.0"
    );
    anyhow::ensure!(args.threads > 0, "threads must be >=1");
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
        constant,
        alphabet,
        args.threads,
    )?;
    eprintln!("# analyse variant groups");
    eprintln!(" . raw variant groups: {}", raw_groups.groups.len());

    let sorted_groups = group_sorting::sort_and_deduplicate_groups(raw_groups, alphabet)?;
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

    let (alen, amiss, aconstant) = output::write_outputs(
        &args.output,
        &res.species_names,
        res.concat,
        res.partitions,
        res.partition_names,
        args.nj,
        args.threads,
    )?;

    eprintln!("# output files");
    eprintln!(
        " . alignment: length={} constant={:.1}% missing={:.1}%",
        alen, aconstant, amiss
    );
    if args.nj {
        eprintln!(" . NJ tree");
    }
    drop(genomes_tmpdir);
    Ok(())
}

fn resolve_alphabet_parameters(
    recode: Option<RecodeScheme>,
    requested_k: Option<usize>,
    requested_constant: Option<usize>,
) -> anyhow::Result<(recode::Alphabet, usize, usize)> {
    let alphabet = recode::Alphabet::new(recode);
    let (default_k, default_constant) = if recode.is_some() { (13, 3) } else { (8, 1) };
    let k = requested_k.unwrap_or(default_k);
    let constant = requested_constant.unwrap_or(default_constant.min(k));
    anyhow::ensure!((1..=alphabet.max_packed_k()).contains(&k), "invalid k");
    anyhow::ensure!(constant <= k, "constant <= k");
    Ok((alphabet, k, constant))
}

#[cfg(test)]
mod argument_tests {
    use super::*;

    #[test]
    fn alphabet_specific_defaults_and_explicit_overrides() {
        let (_, k, c) = resolve_alphabet_parameters(None, None, None).unwrap();
        assert_eq!((k, c), (8, 1));
        let (_, k, c) = resolve_alphabet_parameters(Some(RecodeScheme::SR6), None, None).unwrap();
        assert_eq!((k, c), (13, 3));
        let (_, k, c) = resolve_alphabet_parameters(None, Some(4), Some(2)).unwrap();
        assert_eq!((k, c), (4, 2));
        let (_, k, c) =
            resolve_alphabet_parameters(Some(RecodeScheme::Dayhoff6), Some(2), None).unwrap();
        assert_eq!((k, c), (2, 2));
    }

    #[test]
    fn k_bounds_follow_the_packed_alphabet_width() {
        assert!(resolve_alphabet_parameters(None, Some(12), None).is_ok());
        assert!(resolve_alphabet_parameters(None, Some(13), None).is_err());
        assert!(resolve_alphabet_parameters(Some(RecodeScheme::KGB6), Some(21), None).is_ok());
        assert!(resolve_alphabet_parameters(Some(RecodeScheme::KGB6), Some(22), None).is_err());
        assert!(resolve_alphabet_parameters(None, Some(2), Some(3)).is_err());
    }
}
