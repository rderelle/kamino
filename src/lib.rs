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
//!
//! Main parameters:
//!   -i, --input-directory <INPUT>        Directory containing proteome FASTA files
//!   -I, --input-file <INPUT_FILE>        TSV table with `species_name<TAB>path_to_fasta` rows
//!   -o, --output <OUTPUT>                Prefix for output files [default: kamino]
//!   -k, --k <K>                          k-mer size used for anchor extraction [k=8]
//!   -f, --min-freq <MIN_FREQ>            Minimum fraction of species present per alignment position [f=0.85]
//!   -c, --constant <CONSTANT>            Number of 'constant' positions added to each partition [c=3]
//!   -l, --length-middle <LENGTH_MIDDLE>  Maximum amino-acid length between two adjacent shared anchors [l=35]
//!   -m, --mask <MASK>                    Consecutive amino-acid differences required for masking [m=5]; 0 disables
//!   -t, --threads <THREADS>              Number of threads [t=1]
//!
//! Optional input:
//!       --genomes  Treat inputs as bacterial genomes and predict proteins first
//!
//! Optional output:
//!       --nj                     Build a neighbor-joining tree
//!   -b, --bootstrap <BOOTSTRAP>  Number of bootstrap replicates; requires --nj
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
//! amino-acid k-mer in each variant group, next to the middle positions. With the
//! default value of c = 3, constant positions represent about 50% of the alignment.
//!
//! The --mask parameter controls the amino-acid masking performed by kamino to
//! prevent long runs of polymorphism from being retained in the final alignment.
//! These runs correspond to genuine but unwanted polymorphisms, such as
//! micro-inversions, or, less frequently, errors made by kamino, such as misaligned
//! paths caused by two consecutive indels. The minimum length of polymorphic runs to
//! be masked can be decreased with this parameter to make the filtering more
//! stringent.
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
//! When `-b/--bootstrap` is supplied together with `--nj`, alignment columns are resampled
//! with replacement and bootstrap percentages are written as internal-node labels in the NJ tree.
//! Please note that bootstrap analyses on large datasets (e.g. >1,000 samples) can be
//! computationally intensive.
//!
use anyhow::Context;
use clap::Parser;

const DEFAULT_K: usize = 8;

mod amino_acid;
mod group_extraction;
mod group_filtering;
mod group_sorting;
mod io;
mod output;
mod phylo;
mod proba_filter;
mod protein_prediction;

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

    /// k-mer size used for anchor extraction [k=8].
    #[arg(short, long, help_heading = "Main parameters")]
    pub k: Option<usize>,

    /// Minimum fraction of species present per alignment position [f=0.85].
    #[arg(
        short = 'f',
        long = "min-freq",
        default_value_t = 0.85,
        hide_default_value = true,
        help_heading = "Main parameters"
    )]
    pub min_freq: f32,

    /// Number of 'constant' positions added to each partition [c=3].
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

    /// Number of bootstrap replicates; requires --nj.
    #[arg(
        short = 'b',
        long = "bootstrap",
        requires = "nj",
        help_heading = "Optional output"
    )]
    pub bootstrap: Option<usize>,
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
    if args.nj {
        parameters.push("nj=true".to_string());
    }
    if let Some(bootstrap) = args.bootstrap {
        parameters.push(format!("bootstrap={bootstrap}"));
    }

    eprintln!("kamino {}", env!("CARGO_PKG_VERSION"));
    eprintln!("parameters: {}", parameters.join(" "));
}

/// Validate arguments, collect inputs, run the analysis, and write output files.
pub fn run_with_args(args: Args) -> anyhow::Result<()> {
    // Defaults and bounds are kept here so tests and the CLI share identical behavior.
    let k = args.k.unwrap_or(DEFAULT_K);
    let constant = args.constant.unwrap_or(3usize.min(k));
    print_startup_banner(&args, k, constant);
    anyhow::ensure!(
        (0.6..=1.0).contains(&args.min_freq),
        "min_freq must be between 0.6 and 1.0"
    );
    anyhow::ensure!(args.threads > 0, "threads must be >=1");
    if let Some(bootstrap) = args.bootstrap {
        anyhow::ensure!(args.nj, "--bootstrap requires --nj");
        anyhow::ensure!(bootstrap > 0, "bootstrap replicates must be >=1");
    }
    anyhow::ensure!(
        (1..=amino_acid::MAX_PACKED_K).contains(&k),
        "invalid k (expected 1..={})",
        amino_acid::MAX_PACKED_K
    );
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
        constant,
        args.threads,
    )?;
    eprintln!("# analyse variant groups");
    eprintln!(" . raw variant groups: {}", raw_groups.groups.len());

    let sorted_groups = group_sorting::sort_and_deduplicate_groups(raw_groups)?;
    eprintln!(
        " . sorted variant groups: {}",
        sorted_groups.raw_candidates.len()
    );

    let res =
        group_filtering::filter_groups(sorted_groups, args.min_freq, args.mask, args.threads)?;
    eprintln!(" . filtered variant groups: {}", res.partitions.len());

    let (alen, amiss) = output::write_outputs(
        &args.output,
        &res.species_names,
        res.concat,
        res.partitions,
        res.partition_names,
        args.nj,
        args.bootstrap,
        args.threads,
    )?;

    eprintln!("# output files");
    eprintln!(" . alignment: length={} missing={:.1}%", alen, amiss);
    if args.nj {
        if let Some(bootstrap) = args.bootstrap {
            eprintln!(" . NJ tree + {bootstrap} bootstrap replicates");
        } else {
            eprintln!(" . NJ tree");
        }
    }
    drop(genomes_tmpdir);
    Ok(())
}

#[cfg(test)]
mod cli_tests {
    use super::*;

    #[test]
    fn default_k_is_selected_by_the_pipeline() {
        let args = Args::try_parse_from(["kamino", "-i", "input"]).unwrap();
        assert_eq!(args.k, None);
        assert_eq!(DEFAULT_K, 8);
    }

    #[test]
    fn k_twelve_parses_and_removed_alphabet_options_are_unknown() {
        assert_eq!(
            Args::try_parse_from(["kamino", "-i", "input", "-k", "12"])
                .unwrap()
                .k,
            Some(12)
        );
        assert!((1..=amino_acid::MAX_PACKED_K).contains(&12));
        assert!(!(1..=amino_acid::MAX_PACKED_K).contains(&13));
        assert!(Args::try_parse_from(["kamino", "-i", "input", "-r", "sr6"]).is_err());
        assert!(Args::try_parse_from(["kamino", "-i", "input", "--recode", "sr6"]).is_err());
    }
}
