//! Output writers for the concatenated alignment, missing-data report, partitions,
//! and optional neighbor-joining tree.
use anyhow::{Context, Result};
use std::ffi::OsStr;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

pub fn output_paths(out_base: &Path) -> (PathBuf, PathBuf, PathBuf, PathBuf) {
    // Derive all output filenames from one prefix so `-o results/foo` writes
    // `results/foo_alignment.fas`, `results/foo_missing.tsv`, and so on.
    fn build_path(base: &Path, suffix: &str, ext: &str) -> PathBuf {
        let stem = base
            .file_stem()
            .or_else(|| base.file_name())
            .unwrap_or_else(|| OsStr::new("output"));
        let mut name = stem.to_os_string();
        name.push(suffix);
        let mut path = base.parent().map(Path::to_path_buf).unwrap_or_default();
        path.push(name);
        path.set_extension(ext);
        path
    }
    (
        build_path(out_base, "_alignment", "fas"),
        build_path(out_base, "_missing", "tsv"),
        build_path(out_base, "_partitions", "tsv"),
        build_path(out_base, "_NJ", "tree"),
    )
}

pub fn write_outputs(
    out_base: &Path,
    species: &[String],
    concat: Vec<Vec<u8>>,
    partitions: Vec<(usize, usize, usize)>,
    partition_names: Vec<String>,
    generate_nj: bool,
    threads: usize,
) -> Result<(usize, f64)> {
    // Calculate overall missingness once before writing the FASTA rows.
    let (fas_path, tsv_path, partitions_path, tree_path) = output_paths(out_base);
    let mut w = BufWriter::new(File::create(&fas_path)?);
    let total_len = concat.first().map(|v| v.len()).unwrap_or(0);
    let total_pos = total_len * species.len();
    let total_missing: usize = concat
        .iter()
        .map(|s| s.iter().filter(|&&b| b == b'-' || b == b'X').count())
        .sum();
    let miss = if total_pos == 0 {
        0.0
    } else {
        100.0 * total_missing as f64 / total_pos as f64
    };
    // FASTA output is wrapped at 60 characters for broad tool compatibility.
    for (sid, name) in species.iter().enumerate() {
        writeln!(w, ">{}", name)?;
        for ch in concat[sid].chunks(60) {
            w.write_all(ch)?;
            w.write_all(b"\n")?;
        }
    }
    w.flush()?;
    // Per-species missingness helps users identify problematic samples.
    let mut tw = BufWriter::new(File::create(&tsv_path)?);
    writeln!(tw, "species\tpct_missing")?;
    for (sid, name) in species.iter().enumerate() {
        let len = concat[sid].len().max(1);
        let m = concat[sid]
            .iter()
            .filter(|&&b| b == b'-' || b == b'X')
            .count();
        writeln!(tw, "{}\t{:.1}", name, 100.0 * m as f64 / len as f64)?;
    }
    tw.flush()?;
    // Partitions use the coordinates returned by the extraction pipeline and keep
    // the chosen consensus protein label next to each block.
    let mut pw = BufWriter::new(
        File::create(&partitions_path).with_context(|| format!("create {:?}", partitions_path))?,
    );
    writeln!(pw, "start_pos\tend_pos\tlength\tconsensus protein name")?;
    for ((s, e, l), n) in partitions.into_iter().zip(partition_names) {
        writeln!(pw, "{}\t{}\t{}\t{}", s, e, l, n)?;
    }
    pw.flush()?;
    if generate_nj {
        // Tree generation is optional because it adds an O(n²·L) distance pass.
        let tree = crate::phylo::nj_tree_newick(species, &concat, threads)
            .map_err(|e| anyhow::anyhow!(e))?;
        let mut tw = BufWriter::new(File::create(&tree_path)?);
        writeln!(tw, "{}", tree)?;
        tw.flush()?;
    }
    Ok((total_len, miss))
}
