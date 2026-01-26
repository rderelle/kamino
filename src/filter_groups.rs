use anyhow::{bail, Result};
use hashbrown::{HashMap, HashSet};

use crate::io::SpeciesInput;
use crate::revert_aminoacid::{self, GroupLayout, KmerNameStats, WordInterner};
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

pub(crate) struct GroupKeptMatrix {
    pub k: usize,
    pub data: Vec<u8>,
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

fn mask_middle_if_diff_run_kept(
    matrix: &mut GroupKeptMatrix,
    layout: &GroupLayout,
    m: usize,
    n: usize,
) {
    let middle_len = layout.middle_len as usize;
    if middle_len == 0 || middle_len < m {
        return;
    }

    let mut consensus: Vec<u8> = vec![b'X'; middle_len];
    for (col, consensus_cell) in consensus.iter_mut().enumerate().take(middle_len) {
        let mut counts = [0usize; 256];
        let mut valid_votes = 0usize;
        for sid in 0..n {
            let idx = sid * matrix.k + col;
            let b = matrix.data[idx];
            if b == b'-' || b == b'X' {
                continue;
            }
            counts[b as usize] += 1;
            valid_votes += 1;
        }
        if valid_votes == 0 {
            *consensus_cell = b'X';
            continue;
        }
        let mut best = b'X';
        let mut best_count = 0usize;
        for (aa, &count) in counts.iter().enumerate() {
            if count > best_count {
                best_count = count;
                best = aa as u8;
            }
        }
        if best_count * 2 >= valid_votes {
            *consensus_cell = best;
        } else {
            *consensus_cell = b'X';
        }
    }

    for sid in 0..n {
        let row_start = sid * matrix.k;
        let row = &mut matrix.data[row_start..row_start + matrix.k];
        let mut run = 0usize;
        let mut should_mask = false;
        for col in 0..middle_len {
            let a = row[col];
            let c = consensus[col];
            if a == b'-' || a == b'X' || c == b'-' || c == b'X' {
                run = 0;
                continue;
            }
            if a != c {
                run += 1;
            } else {
                run = 0;
            }
            if run >= m {
                should_mask = true;
                break;
            }
        }
        if should_mask {
            for cell in row.iter_mut().take(middle_len) {
                *cell = b'X';
            }
        }
    }
}

fn filter_group_kept(
    matrix: &mut GroupKeptMatrix,
    layout: &GroupLayout,
    bubble_ratio: f32,
    max_middle_len: usize,
    mask_m: usize,
    n: usize,
) -> FilterOutcome {
    if matrix.k == 0 || matrix.data.is_empty() {
        return FilterOutcome::Skipped;
    }

    mask_middle_if_diff_run_kept(matrix, layout, mask_m, n);

    let middle_len = layout.middle_len as usize;
    let tail_keep = layout.tail_const_keep as usize;
    let missing_cutoff = 1.0 - bubble_ratio as f64;

    let mut keep_tail_cols: Vec<usize> = Vec::with_capacity(tail_keep);
    for col in 0..tail_keep {
        let mut missing = 0usize;
        let idx = middle_len + col;
        for sid in 0..n {
            let b = matrix.data[sid * matrix.k + idx];
            if b == b'-' || b == b'X' {
                missing += 1;
            }
        }
        let missing_ratio = (missing as f64) / (n as f64);
        if missing_ratio <= missing_cutoff {
            keep_tail_cols.push(idx);
        }
    }

    let mut keep_middle_cols: Vec<usize> = Vec::with_capacity(middle_len);
    let mut has_polymorphic_middle = false;

    for col in 0..middle_len {
        let mut counts = [0usize; 256];
        let mut missing = 0usize;

        for sid in 0..n {
            let b = matrix.data[sid * matrix.k + col];
            if b == b'-' || b == b'X' {
                missing += 1;
            }
            counts[b as usize] += 1;
        }

        let missing_ratio = (missing as f64) / (n as f64);
        if missing_ratio > missing_cutoff {
            continue;
        }

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

        keep_middle_cols.push(col);
    }

    if keep_middle_cols.is_empty() {
        return FilterOutcome::DroppedMiddle;
    }

    if !has_polymorphic_middle {
        return FilterOutcome::DroppedMiddle;
    }

    if keep_middle_cols.len() > max_middle_len {
        return FilterOutcome::DroppedLength;
    }

    let block_len_filtered = keep_middle_cols.len() + keep_tail_cols.len();
    debug_assert!(block_len_filtered > 0);

    let mut filtered_rows: Vec<Vec<u8>> = Vec::with_capacity(n);
    for sid in 0..n {
        let row_start = sid * matrix.k;
        let row = &matrix.data[row_start..row_start + matrix.k];
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
    inputs: &[SpeciesInput],
    species: &[String],
    recode_scheme: crate::RecodeScheme,
    groups: &VariantGroups,
    layouts: &[GroupLayout],
    bubble_ratio: f32,
    max_middle_len: usize,
    mask_m: usize,
    scan_k: usize,
    num_threads: usize,
) -> Result<VariantGroupAlignment> {
    let _ = num_threads;
    let n = species.len();

    // Per-species buffers for the final concatenated alignment.
    let mut concat: Vec<Vec<u8>> = vec![Vec::new(); n];
    let mut partitions: Vec<(usize, usize, usize)> = Vec::new();
    let mut partition_names: Vec<String> = Vec::new();
    let mut current_pos: usize = 0;

    let mut dropped_blocks_due_to_middle_filter = 0usize;
    let mut dropped_blocks_due_to_length = 0usize;

    if inputs.len() != species.len() {
        bail!(
            "Input FASTA count ({}) does not match graph species count ({}).",
            inputs.len(),
            species.len()
        );
    }

    let mut name_to_sid: HashMap<String, usize> = HashMap::new();
    for (sid, name) in species.iter().enumerate() {
        name_to_sid.insert(name.clone(), sid);
    }

    let mut kmers_of_interest: HashSet<u64> = HashSet::new();
    for layout in layouts {
        for &km in &layout.unique_kmers {
            kmers_of_interest.insert(km);
        }
    }

    let mut matrices: Vec<GroupKeptMatrix> = layouts
        .iter()
        .map(|layout| GroupKeptMatrix {
            k: layout.kept_len(),
            data: vec![b'-'; layout.kept_len() * n],
        })
        .collect();

    let mut group_name_stats: Vec<KmerNameStats> = vec![KmerNameStats::default(); layouts.len()];
    let mut word_interner = WordInterner::new();

    for input in inputs {
        let Some(&sid) = name_to_sid.get(&input.name) else {
            bail!(
                "Species {:?} from input list not found in graph.",
                input.name
            );
        };

        let species_map = revert_aminoacid::build_species_kmer_map(
            input,
            recode_scheme,
            scan_k,
            &kmers_of_interest,
            &mut word_interner,
        )?;

        let mut seen_kmers: HashSet<u64> = HashSet::new();

        for (gidx, layout) in layouts.iter().enumerate() {
            if layout.paths.is_empty() || layout.kept_len() == 0 {
                continue;
            }

            let matrix = &mut matrices[gidx];
            let row_start = sid * matrix.k;
            let row = &mut matrix.data[row_start..row_start + matrix.k];

            seen_kmers.clear();

            let paths = &groups[&layout.key];
            for (path_layout, &path_idx) in layout.paths.iter().zip(layout.path_indices.iter()) {
                let path = &paths[path_idx];
                if !path.species.get(sid).map(|b| *b).unwrap_or(false) {
                    continue;
                }

                for (kmer, pos0) in path_layout.kmers.iter().zip(path_layout.pos0.iter()) {
                    let Some(&idx) = species_map.map.get(kmer) else {
                        continue;
                    };

                    let offset = idx as usize * scan_k;
                    let window = &species_map.windows[offset..offset + scan_k];

                    for (t, &newc) in window.iter().enumerate().take(scan_k) {
                        let pos = *pos0 as usize + t;
                        if pos >= layout.bl as usize {
                            break;
                        }
                        let Some(kept_idx) = layout.kept_index(pos as u32) else {
                            continue;
                        };

                        let cur = row[kept_idx];
                        if cur == b'-' {
                            row[kept_idx] = newc;
                        } else if newc == b'-' {
                            // keep cur
                        } else if cur == newc {
                            // identical call; keep cur
                        } else {
                            row[kept_idx] = b'X';
                        }
                    }

                    if seen_kmers.insert(*kmer) {
                        let name_stats = &species_map.name_stats[idx as usize];
                        if name_stats.total_sets > 0 {
                            let group_stats = &mut group_name_stats[gidx];
                            group_stats.total_sets =
                                group_stats.total_sets.saturating_add(name_stats.total_sets);
                            for (word_id, count) in &name_stats.word_counts {
                                if let Some(entry) = group_stats
                                    .word_counts
                                    .iter_mut()
                                    .find(|(existing_id, _)| existing_id == word_id)
                                {
                                    entry.1 = entry.1.saturating_add(*count);
                                } else {
                                    group_stats.word_counts.push((*word_id, *count));
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    let word_list = word_interner.words;

    for (gidx, (layout, matrix)) in layouts.iter().zip(matrices.iter_mut()).enumerate() {
        match filter_group_kept(matrix, layout, bubble_ratio, max_middle_len, mask_m, n) {
            FilterOutcome::Kept {
                filtered_rows,
                filtered_len,
            } => {
                let start_pos = current_pos;
                let end_pos = start_pos + filtered_len - 1;
                partitions.push((start_pos, end_pos, filtered_len));
                partition_names.push(build_consensus_name(&group_name_stats[gidx], &word_list));
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

    Ok(VariantGroupAlignment {
        concat,
        partitions,
        partition_names,
        dropped_blocks_due_to_middle_filter,
        dropped_blocks_due_to_length,
    })
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
