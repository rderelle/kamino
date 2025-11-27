use anyhow::{Context, Result};
use std::ffi::OsStr;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::graph::Graph;
use crate::traverse::VariantGroups;
use crate::variant_groups::{
    build_group_block, build_variant_group_alignment, GroupOutcome, VariantGroupAlignment,
};
use hashbrown::HashMap;
use rayon::prelude::*;

/// Assemble the final concatenated alignment matrix for all variant groups.
///
/// The function walks each variant group in deterministic order, applies the
/// middle-filter rules, and stitches the retained blocks into a single matrix
/// while tracking the block boundaries. Consensus repainting relies on the
/// per-species k-mer consensus gathered earlier in the pipeline.
#[allow(clippy::too_many_arguments)]
pub(crate) fn build_concatenated_alignment(
    groups: &VariantGroups,
    k: usize,
    k1: usize,
    head_keep: usize,
    sym_bits: u32,
    scan_k: usize,
    species_kmer_maps: &[HashMap<u64, usize>],
    species_kmer_consensus: &[Vec<Option<Vec<u8>>>],
    bubble_ratio: f32,
    max_middle_len: usize,
    n: usize,
) -> Result<VariantGroupAlignment> {
    // Per-species buffers for the final concatenated alignment.
    let mut concat: Vec<Vec<u8>> = vec![Vec::new(); n];
    let mut partitions: Vec<(usize, usize, usize)> = Vec::new();
    let mut current_pos: usize = 0;
    let bubble_ratio_f = bubble_ratio as f64;

    // Deterministic iteration: sort (start,end)
    let mut keys: Vec<(u64, u64)> = groups.keys().copied().collect();
    keys.sort_unstable();

    let mut dropped_blocks_due_to_middle_filter = 0usize;
    let mut dropped_blocks_due_to_length = 0usize;

    let group_results: Vec<GroupOutcome> = keys
        .par_iter()
        .map(|&key| {
            let paths = &groups[&key];
            build_group_block(
                key,
                paths,
                k,
                k1,
                head_keep,
                sym_bits,
                scan_k,
                species_kmer_maps,
                species_kmer_consensus,
                bubble_ratio_f,
                max_middle_len,
                n,
            )
        })
        .collect();

    for outcome in group_results {
        match outcome {
            GroupOutcome::Kept {
                filtered_rows,
                filtered_len,
            } => {
                let start_pos = current_pos;
                let end_pos = start_pos + filtered_len - 1;
                partitions.push((start_pos, end_pos, filtered_len));
                current_pos += filtered_len;

                for (sid, row) in filtered_rows.into_iter().enumerate() {
                    let dst = &mut concat[sid];
                    dst.reserve(row.len());
                    dst.extend_from_slice(&row);
                }
            }
            GroupOutcome::DroppedMiddle => {
                dropped_blocks_due_to_middle_filter += 1;
            }
            GroupOutcome::DroppedLength => {
                dropped_blocks_due_to_length += 1;
            }
            GroupOutcome::Skipped => {}
        }
    }

    Ok(VariantGroupAlignment {
        concat,
        partitions,
        dropped_blocks_due_to_middle_filter,
        dropped_blocks_due_to_length,
    })
}

pub fn output_paths(out_base: &Path) -> (PathBuf, PathBuf, PathBuf) {
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
    (fas, tsv, partitions)
}

/// Write FASTA and TSV outputs while trimming the head and middle segments.
///
/// The function orchestrates variant-group alignment building, applies user
/// limits on leading positions and middle length, and emits the three output
/// files alongside summary stats.
#[allow(clippy::too_many_arguments)]
pub fn write_outputs_with_head(
    input_dir: &Path,
    out_base: &Path,
    g: &Graph,
    groups: &VariantGroups,
    k: usize,
    head_max: usize,   // user-defined n; must be ≤ k-1
    bubble_ratio: f32, // proportion of species present required to keep a column
    max_middle_len: usize,
    num_threads: usize,
) -> Result<(usize, f64)> {
    let species = &g.species_names;

    let alignment = build_variant_group_alignment(
        input_dir,
        g,
        groups,
        k,
        head_max,
        bubble_ratio,
        max_middle_len,
        num_threads,
    )?;

    let (fas_path, tsv_path, partitions_path) = output_paths(out_base);
    let (concat, partitions, dropped_middle, dropped_length) = (
        alignment.concat,
        alignment.partitions,
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
    writeln!(pw, "start_pos\tend_pos\tlength")?;
    for (start, end, len) in partitions {
        writeln!(pw, "{}\t{}\t{}", start, end, len)?;
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
    Ok((total_len, alignment_missing_pct))
}
