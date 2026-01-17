use hashbrown::HashMap;

use crate::middle_mask::mask_middle_if_diff_run;
use crate::recode::RECODE_BITS_PER_SYMBOL;
use crate::revert_aminoacid::{self, KmerNameStats};
use crate::traverse::VariantGroups;

/// Alignment output pieces produced for all variant groups.
///
/// * `concat` stores the stitched alignment rows, one per species.
/// * `partitions` records the start/end indices for each variant group block.
/// * `dropped_blocks_due_to_middle_filter` counts blocks discarded by the
///   uniqueness/missing filter on middle positions.
/// * `dropped_blocks_due_to_length` counts blocks that were too long once the
///   middle segment limit was applied.
pub struct VariantGroupAlignment {
    pub concat: Vec<Vec<u8>>,
    pub partitions: Vec<(usize, usize, usize)>,
    pub partition_names: Vec<String>,
    pub dropped_blocks_due_to_middle_filter: usize,
    pub dropped_blocks_due_to_length: usize,
}

pub(crate) enum FilterOutcome {
    Kept {
        filtered_rows: Vec<Vec<u8>>,
        filtered_len: usize,
    },
    DroppedMiddle,
    DroppedLength,
    Skipped,
}

pub(crate) fn filter_group_block(
    block: &mut [Vec<u8>],
    k1: usize,
    head_keep: usize,
    bubble_ratio: f32,
    max_middle_len: usize,
    mask_m: usize,
) -> FilterOutcome {
    let n = block.len();
    let Some(bl) = block.first().map(|row| row.len()) else {
        return FilterOutcome::Skipped;
    };
    if bl == 0 {
        return FilterOutcome::Skipped;
    }

    // Region indices
    let lead_full = k1;
    let trail_full = k1;

    let middle_start = lead_full;
    let middle_end = bl - trail_full; // exclusive
    let middle_len = middle_end - middle_start;

    mask_middle_if_diff_run(block, k1, mask_m);

    // Missing cutoff based on presence ratio requirement (presence_ratio = 1 - missing_ratio >= bubble_ratio)
    let missing_cutoff = 1.0 - bubble_ratio as f64;

    // Trailing constant slice: the end k-mer contributes k-1 columns
    let tail_const_start = bl - k1;

    let tail_const_keep = head_keep.min(bl - tail_const_start);
    let tail_const_end = tail_const_start.saturating_add(tail_const_keep);

    // Columns from the trailing context (beginning of end k-mer determined from
    // the last non-variable column): filter ONLY on missing ratio (constant
    // columns allowed).
    let mut keep_tail_cols: Vec<usize> = Vec::with_capacity(tail_const_keep);
    for col in tail_const_start..tail_const_end {
        let mut missing = 0usize;
        for row in block.iter().take(n) {
            let b = row[col];
            if b == b'-' || b == b'X' {
                missing += 1;
            }
        }
        let missing_ratio = (missing as f64) / (n as f64);
        if missing_ratio <= missing_cutoff {
            keep_tail_cols.push(col);
        }
    }

    // Filter middle positions on missing ratio (polymorphism is checked at the group level)
    let mut keep_middle_cols: Vec<usize> = Vec::with_capacity(middle_len);
    let mut has_polymorphic_middle = false;

    for col in middle_start..middle_end {
        if col >= tail_const_start && col < tail_const_end {
            // Skip positions reserved for the trailing constant context.
            continue;
        }

        let mut counts = [0usize; 256];
        let mut missing = 0usize;

        for row in block.iter().take(n) {
            let b = row[col];
            if b == b'-' || b == b'X' {
                missing += 1;
            }
            counts[b as usize] += 1;
        }

        let missing_ratio = (missing as f64) / (n as f64);
        if missing_ratio > missing_cutoff {
            // Too much missing: drop this column.
            continue;
        }

        // Check if this column is polymorphic among non-missing AAs.
        let mut unique_non_gap = 0usize;
        for (aa, &cnt) in counts.iter().enumerate() {
            if cnt == 0 {
                continue;
            }
            if aa != b'-' as usize && aa != b'X' as usize {
                unique_non_gap += 1;
                if unique_non_gap > 1 {
                    has_polymorphic_middle = true;
                    break;
                }
            }
        }

        // Keep all columns that pass the missing filter, even if monomorphic.
        keep_middle_cols.push(col);
    }

    // If nothing survives the missing filter, drop the block.
    if keep_middle_cols.is_empty() {
        return FilterOutcome::DroppedMiddle;
    }

    // If no middle column is polymorphic, drop the block.
    if !has_polymorphic_middle {
        return FilterOutcome::DroppedMiddle;
    }

    if keep_middle_cols.len() > max_middle_len {
        return FilterOutcome::DroppedLength;
    }

    let block_len_filtered = keep_middle_cols.len() + keep_tail_cols.len();
    debug_assert!(block_len_filtered > 0);

    let mut filtered_rows: Vec<Vec<u8>> = Vec::with_capacity(n);
    for row in block.iter().take(n) {
        let mut dst: Vec<u8> = Vec::with_capacity(block_len_filtered);

        for &col in &keep_middle_cols {
            dst.push(row[col]);
        }
        for &col in &keep_tail_cols {
            dst.push(row[col]);
        }
        filtered_rows.push(dst);
    }

    FilterOutcome::Kept {
        filtered_rows,
        filtered_len: block_len_filtered,
    }
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn build_concatenated_alignment_streaming(
    groups: &VariantGroups,
    k: usize,
    head_keep: usize,
    bubble_ratio: f32,
    max_middle_len: usize,
    mask_m: usize,
    scan_k: usize,
    species_kmer_maps: &[HashMap<u64, usize>],
    species_kmer_consensus: &[Vec<Option<Vec<u8>>>],
    species_kmer_name_stats: &[Vec<KmerNameStats>],
    word_list: &[String],
    n: usize,
) -> VariantGroupAlignment {
    let k1 = k - 1;
    let sym_bits = RECODE_BITS_PER_SYMBOL;

    // Per-species buffers for the final concatenated alignment.
    let mut concat: Vec<Vec<u8>> = vec![Vec::new(); n];
    let mut partitions: Vec<(usize, usize, usize)> = Vec::new();
    let mut partition_names: Vec<String> = Vec::new();
    let mut current_pos: usize = 0;

    // Deterministic iteration: sort (start,end)
    let mut keys: Vec<(u64, u64)> = groups.keys().copied().collect();
    keys.sort_unstable();

    let mut dropped_blocks_due_to_middle_filter = 0usize;
    let mut dropped_blocks_due_to_length = 0usize;

    for key in keys {
        let paths = &groups[&key];
        let Some((mut block, group_name_stats)) = revert_aminoacid::build_raw_group_block(
            key,
            paths,
            k,
            k1,
            head_keep,
            sym_bits,
            scan_k,
            species_kmer_maps,
            species_kmer_consensus,
            species_kmer_name_stats,
            n,
        ) else {
            continue;
        };

        match filter_group_block(
            &mut block,
            k1,
            head_keep,
            bubble_ratio,
            max_middle_len,
            mask_m,
        ) {
            FilterOutcome::Kept {
                filtered_rows,
                filtered_len,
            } => {
                let start_pos = current_pos;
                let end_pos = start_pos + filtered_len - 1;
                partitions.push((start_pos, end_pos, filtered_len));
                partition_names.push(build_consensus_name(&group_name_stats, word_list));
                current_pos += filtered_len;

                for (sid, row) in filtered_rows.into_iter().enumerate() {
                    let dst = &mut concat[sid];
                    dst.reserve(row.len());
                    dst.extend_from_slice(&row);
                }
            }
            FilterOutcome::DroppedMiddle => {
                dropped_blocks_due_to_middle_filter += 1;
            }
            FilterOutcome::DroppedLength => {
                dropped_blocks_due_to_length += 1;
            }
            FilterOutcome::Skipped => {}
        }
    }

    VariantGroupAlignment {
        concat,
        partitions,
        partition_names,
        dropped_blocks_due_to_middle_filter,
        dropped_blocks_due_to_length,
    }
}

fn build_consensus_name(stats: &KmerNameStats, word_list: &[String]) -> String {
    if stats.total_sets == 0 {
        return String::new();
    }
    let mut word_ids: Vec<u32> = stats
        .word_counts
        .iter()
        .filter(|(_, count)| (*count as u64) * 2 >= stats.total_sets as u64)
        .map(|(word_id, _)| *word_id)
        .collect();
    word_ids.sort_unstable();
    let mut output = String::new();
    for (idx, word_id) in word_ids.iter().enumerate() {
        let Some(word) = word_list.get(*word_id as usize) else {
            continue;
        };
        if idx > 0 && !output.is_empty() {
            output.push(' ');
        }
        output.push_str(word);
    }
    output
}
