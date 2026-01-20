use anyhow::{Context, Result};
use std::ffi::OsStr;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::filter_groups::build_concatenated_alignment_streaming;
use crate::graph::Graph;
use crate::io::SpeciesInput;
use crate::phylo;
use crate::revert_aminoacid;
use crate::traverse::VariantGroups;

pub fn output_paths(out_base: &Path) -> (PathBuf, PathBuf, PathBuf, PathBuf) {
    /// Build an output path next to the base path while keeping the stem.
    fn build_path(base: &Path, suffix: &str, ext: &str) -> PathBuf {
        let stem = base
            .file_stem()
            .or_else(|| base.file_name())
            .unwrap_or_else(|| OsStr::new("output"));
        let mut name = stem.to_os_string();
        name.push(suffix);

        let mut path = if let Some(parent) = base.parent() {
            let p = parent.to_path_buf();
            if !p.as_os_str().is_empty() {
                p
            } else {
                PathBuf::new()
            }
        } else {
            PathBuf::new()
        };
        path.push(name);
        path.set_extension(ext);
        path
    }

    let fas = build_path(out_base, "_alignment", "fas");
    let tsv = build_path(out_base, "_missing", "tsv");
    let partitions = build_path(out_base, "_partitions", "tsv");
    let tree = build_path(out_base, "_NJ", "tree");
    (fas, tsv, partitions, tree)
}

/// Write FASTA and TSV outputs while trimming the head and middle segments.
///
/// The function orchestrates variant-group alignment building, applies user
/// limits on leading positions and middle length, and emits the three output
/// files alongside summary stats.
#[allow(clippy::too_many_arguments)]
pub fn write_outputs_with_head(
    inputs: &[SpeciesInput],
    out_base: &Path,
    g: &Graph,
    groups: &VariantGroups,
    k: usize,
    head_max: usize,   // user-defined n; must be ≤ k-1
    bubble_ratio: f32, // proportion of species present required to keep a column
    max_middle_len: usize,
    mask_m: usize,
    num_threads: usize,
    generate_nj: bool,
) -> Result<(usize, f64)> {
    let species = &g.species_names;
    let n = g.n_species;

    let (scan_k, species_kmer_maps, species_kmer_consensus, species_kmer_name_stats, word_list) =
        revert_aminoacid::build_species_kmer_consensus(
            inputs,
            g,
            groups,
            k,
            head_max,
            num_threads,
        )?;

    let alignment = build_concatenated_alignment_streaming(
        groups,
        k,
        head_max,
        bubble_ratio,
        max_middle_len,
        mask_m,
        scan_k,
        &species_kmer_maps,
        &species_kmer_consensus,
        &species_kmer_name_stats,
        &word_list,
        n,
    );

    let (fas_path, tsv_path, partitions_path, tree_path) = output_paths(out_base);
    let (concat, partitions, partition_names, dropped_middle, dropped_length) = (
        alignment.concat,
        alignment.partitions,
        alignment.partition_names,
        alignment.dropped_blocks_due_to_middle_filter,
        alignment.dropped_blocks_due_to_length,
    );

    // FASTA
    let fh = File::create(&fas_path).with_context(|| format!("create {:?}", fas_path))?;
    let mut w = BufWriter::new(fh);

    let total_len: usize = concat.first().map(|v| v.len()).unwrap_or(0);
    let total_positions: usize = total_len.saturating_mul(species.len());
    let total_missing: usize = concat
        .iter()
        .map(|seq| seq.iter().filter(|&&b| b == b'-' || b == b'X').count())
        .sum();
    let alignment_missing_pct = if total_positions == 0 {
        0.0
    } else {
        (total_missing as f64) * 100.0 / (total_positions as f64)
    };

    for (sid, name) in species.iter().enumerate() {
        writeln!(w, ">{}", name)?;
        let row = &concat[sid];
        for chunk in row.chunks(60) {
            w.write_all(chunk)?;
            w.write_all(b"\n")?;
        }
    }
    w.flush()?;

    // TSV with % missing per species
    let th = File::create(&tsv_path).with_context(|| format!("create {:?}", tsv_path))?;
    let mut tw = BufWriter::new(th);
    writeln!(tw, "species\tpct_missing")?;

    for (sid, name) in species.iter().enumerate() {
        let seq = &concat[sid];
        let len = seq.len().max(total_len);
        let missing = seq.iter().filter(|&&b| b == b'-' || b == b'X').count();
        let pct = if len == 0 {
            0.0
        } else {
            (missing as f64) * 100.0 / (len as f64)
        };
        writeln!(tw, "{}\t{:.1}", name, pct)?;
    }
    tw.flush()?;

    // TSV with partition start/end/length information
    let ph =
        File::create(&partitions_path).with_context(|| format!("create {:?}", partitions_path))?;
    let mut pw = BufWriter::new(ph);
    writeln!(pw, "start_pos\tend_pos\tlength\tconsensus protein name")?;
    for ((start, end, len), consensus_name) in
        partitions.into_iter().zip(partition_names.into_iter())
    {
        writeln!(pw, "{}\t{}\t{}\t{}", start, end, len, consensus_name)?;
    }
    pw.flush()?;

    eprintln!(
        "dropped due to unique/missing middle filtering: {}",
        dropped_middle
    );
    eprintln!(
        "dropped due to exceeding middle length limit: {}",
        dropped_length
    );

    // NJ tree
    if generate_nj {
        let tree = phylo::nj_tree_newick(species, &concat, num_threads)
            .map_err(|err| anyhow::anyhow!(err))?;
        let tree_file =
            File::create(&tree_path).with_context(|| format!("create {:?}", tree_path))?;
        let mut tree_writer = BufWriter::new(tree_file);
        writeln!(tree_writer, "{}", tree)?;
        tree_writer.flush()?;
    }

    Ok((total_len, alignment_missing_pct))
}
