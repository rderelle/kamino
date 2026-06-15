//! Final raw group filtering and alignment construction.
use anyhow::Result;
use hashbrown::{HashMap, HashSet};
use rayon::{prelude::*, ThreadPoolBuilder};

use crate::group_extraction::{BlockArena, BlockObs};
use crate::group_sorting::{consensus_rows_by_isolate, RawDirectCandidate, SortedGroups};
use crate::recode::recode_byte;
use crate::RecodeScheme;

/// Fully assembled result returned to the CLI layer before writing files.
pub(crate) struct AlignmentResult {
    /// Species names in the same row order as `concat`.
    pub(crate) species_names: Vec<String>,
    /// Concatenated alignment rows, one byte vector per species.
    pub(crate) concat: Vec<Vec<u8>>,
    /// Partition intervals as `(start, end, length)` in zero-based inclusive coordinates.
    pub(crate) partitions: Vec<(usize, usize, usize)>,
    /// Consensus protein label for each partition.
    pub(crate) partition_names: Vec<String>,
}

struct CandidateGroup {
    /// Flattened per-species amino-acid rows after retained-column projection.
    rows: Vec<u8>,
    /// Number of retained columns in each flattened row.
    kept_len: usize,
    /// Consensus protein name for the partition table.
    name: String,
}

#[inline]
fn is_missing_or_ambiguous(b: u8) -> bool {
    b == b'-' || b == b'X'
}

fn majority_recoded_consensus(
    rows: &[Vec<u8>],
    len: usize,
    recode_scheme: RecodeScheme,
) -> Vec<Option<u8>> {
    let mut consensus = vec![None; len];
    for c in 0..len {
        let mut counts = [0usize; 6];
        for row in rows {
            let b = row[c];
            if is_missing_or_ambiguous(b) {
                continue;
            }
            let recoded = recode_byte(b, recode_scheme);
            if recoded != 255 {
                counts[recoded as usize] += 1;
            }
        }

        let mut best_state = None;
        let mut best_count = 0usize;
        for (state, &count) in counts.iter().enumerate() {
            if count > best_count {
                best_state = Some(state as u8);
                best_count = count;
            }
        }
        consensus[c] = best_state;
    }
    consensus
}

fn mask_rows_with_long_divergence(
    rows: &mut [Vec<u8>],
    consensus: &[Option<u8>],
    mask: usize,
    recode_scheme: RecodeScheme,
) {
    if mask == 0 {
        return;
    }

    for row in rows {
        let mut consecutive_diffs = 0usize;
        let mut should_mask = false;
        for (&row_aa, &consensus_state) in row.iter().zip(consensus) {
            let Some(consensus_state) = consensus_state else {
                consecutive_diffs = 0;
                continue;
            };

            if is_missing_or_ambiguous(row_aa) {
                consecutive_diffs = 0;
                continue;
            }

            let row_state = recode_byte(row_aa, recode_scheme);
            if row_state == 255 {
                consecutive_diffs = 0;
                continue;
            }

            if row_state == consensus_state {
                consecutive_diffs = 0;
            } else {
                consecutive_diffs += 1;
                if consecutive_diffs >= mask {
                    should_mask = true;
                    break;
                }
            }
        }
        if should_mask {
            row.fill(b'X');
        }
    }
}

fn is_polymorphic_column(rows: &[Vec<u8>], col: usize, recode_scheme: RecodeScheme) -> bool {
    let mut first_observed = None;
    for row in rows {
        let b = row[col];
        if is_missing_or_ambiguous(b) {
            continue;
        }
        let recoded = recode_byte(b, recode_scheme);
        if recoded == 255 {
            continue;
        }
        match first_observed {
            None => first_observed = Some(recoded),
            Some(first) if first != recoded => return true,
            Some(_) => {}
        }
    }
    false
}

fn missing_fraction(rows: &[Vec<u8>], col: usize, n_species: usize) -> f32 {
    let miss = rows
        .iter()
        .filter(|r| is_missing_or_ambiguous(r[col]))
        .count();
    (miss as f32) / (n_species as f32)
}

fn consensus_protein_name(observations: &[BlockObs], protein_names: &[String]) -> String {
    // Build a majority-rule label from all non-empty protein names supporting this
    // variant group. A word is retained when it occurs in at least 50% of names.
    // Retained words are then ordered by their average position across the original
    // names in which they occur, so prefixes/suffixes keep their natural location
    // even when the first observed protein name has a different wording.
    let names: Vec<&str> = observations
        .iter()
        .filter_map(|obs| protein_names.get(obs.protein_id as usize))
        .filter_map(|name| {
            name.trim()
                .split_once(char::is_whitespace)
                .map(|(_, desc)| desc.trim())
        })
        .filter(|name| !name.is_empty())
        .collect();

    if names.is_empty() {
        return String::new();
    }

    #[derive(Clone, Copy, Debug)]
    struct WordStats {
        count: usize,
        position_sum: usize,
        first_seen_order: usize,
    }

    let mut word_stats: HashMap<&str, WordStats> = HashMap::new();
    let mut next_seen_order = 0usize;

    for name in &names {
        // Count each word at most once per protein name, matching the majority-rule
        // definition across names rather than across repeated words within one name.
        // If a word appears more than once in the same name, keep its first position.
        let mut seen_in_name: HashSet<&str> = HashSet::new();
        for (position, word) in name.split_whitespace().enumerate() {
            if !seen_in_name.insert(word) {
                continue;
            }

            word_stats
                .entry(word)
                .and_modify(|stats| {
                    stats.count += 1;
                    stats.position_sum += position;
                })
                .or_insert_with(|| {
                    let stats = WordStats {
                        count: 1,
                        position_sum: position,
                        first_seen_order: next_seen_order,
                    };
                    next_seen_order += 1;
                    stats
                });
        }
    }

    let mut majority_words: Vec<(&str, WordStats)> = word_stats
        .into_iter()
        .filter(|(_, stats)| stats.count.saturating_mul(2) >= names.len())
        .collect();

    majority_words.sort_unstable_by(|(word_a, stats_a), (word_b, stats_b)| {
        // Compare average positions exactly as fractions:
        // stats_a.position_sum / stats_a.count vs stats_b.position_sum / stats_b.count.
        (stats_a.position_sum * stats_b.count)
            .cmp(&(stats_b.position_sum * stats_a.count))
            .then_with(|| stats_a.first_seen_order.cmp(&stats_b.first_seen_order))
            .then_with(|| word_a.cmp(word_b))
    });

    majority_words
        .into_iter()
        .map(|(word, _)| word)
        .collect::<Vec<_>>()
        .join(" ")
}

#[allow(clippy::too_many_arguments)]
fn build_candidate_group(
    raw: RawDirectCandidate,
    arena: &BlockArena,
    protein_names: &[String],
    n: usize,
    k: usize,
    constant: usize,
    min_freq: f32,
    mask: usize,
    recode_scheme: RecodeScheme,
) -> Option<CandidateGroup> {
    let name = consensus_protein_name(&raw.selected, protein_names);
    let mut rows = consensus_rows_by_isolate(&raw.selected, arena, n, raw.raw_len);
    let mid_start = k;
    let mid_end = raw.raw_len - k;
    if mid_end <= mid_start {
        return None;
    }
    if mask > 0 {
        let consensus = majority_recoded_consensus(&rows, raw.raw_len, recode_scheme);
        mask_rows_with_long_divergence(&mut rows, &consensus, mask, recode_scheme);
    }

    // Retain only original middle columns satisfying the missing-data
    // frequency threshold, then trim that retained middle to the first and
    // last polymorphic columns in recoded amino-acid space. Constant
    // right-anchor context is appended afterwards and cannot determine or
    // rescue middle polymorphism.
    let mut retained_middle = Vec::new();
    for c in mid_start..mid_end {
        if missing_fraction(&rows, c, n) <= 1.0 - min_freq {
            retained_middle.push(c);
        }
    }
    if retained_middle.is_empty() {
        return None;
    }
    let first_poly_idx = retained_middle
        .iter()
        .position(|&c| is_polymorphic_column(&rows, c, recode_scheme));
    let last_poly_idx = retained_middle
        .iter()
        .rposition(|&c| is_polymorphic_column(&rows, c, recode_scheme));
    let (Some(first_poly_idx), Some(last_poly_idx)) = (first_poly_idx, last_poly_idx) else {
        return None;
    };
    let mut keep_cols: Vec<usize> = retained_middle[first_poly_idx..=last_poly_idx].to_vec();
    
    let tail_keep = constant.min(k);
    for c in mid_end..(mid_end + tail_keep).min(raw.raw_len) {
        if missing_fraction(&rows, c, n) <= 1.0 - min_freq {
            keep_cols.push(c);
        }
    }
    
    let kept_len = keep_cols.len();
    let mut kept_rows = Vec::with_capacity(n * kept_len);
    for row in &rows {
        kept_rows.extend(keep_cols.iter().map(|&c| row[c]));
    }

    Some(CandidateGroup {
        rows: kept_rows,
        kept_len,
        name,
    })
}

pub(crate) fn filter_groups(
    sorted: SortedGroups,
    min_freq: f32,
    constant: usize,
    mask: usize,
    num_threads: usize,
) -> Result<AlignmentResult> {
    let SortedGroups {
        species_names,
        arena,
        raw_candidates,
        protein_names,
        k,
        min_needed: _min_needed,
        n_species: n,
        recode_scheme,
    } = sorted;
    let mut partitions = Vec::new();
    let mut names = Vec::new();
    let mut pos = 0usize;
    // Stage 7: build consensus rows only for non-overlapping raw groups, then
    // apply biological/column filters and construct final alignment candidates.
    // Per-group filtering is independent, but the accepted candidates must be
    // restored to raw-candidate order before the order-dependent concatenation
    // and partition construction below.
    let n_threads = num_threads.max(1);
    let pool = ThreadPoolBuilder::new().num_threads(n_threads).build()?;
    let mut candidates: Vec<(usize, CandidateGroup)> = pool.install(|| {
        raw_candidates
            .into_par_iter()
            .enumerate()
            .filter_map(|(idx, raw)| {
                build_candidate_group(
                    raw,
                    &arena,
                    &protein_names,
                    n,
                    k,
                    constant,
                    min_freq,
                    mask,
                    recode_scheme,
                )
                .map(|candidate| (idx, candidate))
            })
            .collect()
    });
    candidates.sort_unstable_by_key(|(idx, _)| *idx);
    let candidates: Vec<CandidateGroup> = candidates
        .into_iter()
        .map(|(_, candidate)| candidate)
        .collect();
    let final_len: usize = candidates.iter().map(|cand| cand.kept_len).sum();
    let mut concat = vec![Vec::with_capacity(final_len); n];
    // Stage 8: append retained columns to every species row and record partitions.
    for cand in candidates {
        let blen = cand.kept_len;
        for (sid, row_out) in concat.iter_mut().enumerate().take(n) {
            let start = sid * blen;
            row_out.extend_from_slice(&cand.rows[start..start + blen]);
        }
        partitions.push((pos, pos + blen - 1, blen));
        names.push(cand.name);
        pos += blen;
    }
    Ok(AlignmentResult {
        species_names,
        concat,
        partitions,
        partition_names: names,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::group_extraction::BlockArena;
    use crate::group_sorting::RawDirectCandidate;

    fn obs(arena: &mut BlockArena, sid: u16, aa: &[u8]) -> BlockObs {
        BlockObs {
            sid,
            protein_id: sid as u32,
            aa: arena.push(aa).unwrap(),
        }
    }

    fn candidate(rows: &[&[u8]], arena: &mut BlockArena) -> RawDirectCandidate {
        RawDirectCandidate {
            key: (1, 2),
            coverage: rows.len(),
            raw_len: rows[0].len(),
            selected: rows
                .iter()
                .enumerate()
                .map(|(sid, aa)| obs(arena, sid as u16, aa))
                .collect(),
        }
    }

    fn sorted(rows: &[&[u8]]) -> SortedGroups {
        sorted_with_k(rows, 1)
    }

    fn sorted_with_k(rows: &[&[u8]], k: usize) -> SortedGroups {
        let mut arena = BlockArena::new();
        let raw_candidates = vec![candidate(rows, &mut arena)];
        SortedGroups {
            arena,
            species_names: (0..rows.len())
                .map(|sid| format!("species_{sid}"))
                .collect(),
            raw_candidates,
            protein_names: vec!["protein".to_string(); rows.len()],
            k,
            min_needed: 1,
            n_species: rows.len(),
            recode_scheme: RecodeScheme::SR6,
        }
    }

    #[test]
    fn majority_recoded_consensus_ignores_missing_invalid_and_ties_by_state_value() {
        let rows = vec![b"AAX".to_vec(), b"DDB".to_vec(), b"D--".to_vec()];

        assert_eq!(
            majority_recoded_consensus(&rows, 3, RecodeScheme::SR6),
            vec![Some(1), Some(0), None]
        );
    }

    #[test]
    fn same_recoded_state_differences_are_not_masked() {
        let consensus = vec![Some(0); 5];
        let mut rows = vec![b"PPPPP".to_vec()];

        mask_rows_with_long_divergence(&mut rows, &consensus, 3, RecodeScheme::SR6);

        assert_eq!(rows[0], b"PPPPP".to_vec());
    }

    #[test]
    fn fewer_than_mask_consecutive_recoded_differences_are_not_masked() {
        let consensus = vec![Some(0); 5];
        let mut rows = vec![b"ADDAA".to_vec()];

        mask_rows_with_long_divergence(&mut rows, &consensus, 3, RecodeScheme::SR6);

        assert_eq!(rows[0], b"ADDAA".to_vec());
    }

    #[test]
    fn exactly_mask_consecutive_recoded_differences_mask_entire_row() {
        let consensus = vec![Some(0); 5];
        let mut rows = vec![b"ADDDA".to_vec()];

        mask_rows_with_long_divergence(&mut rows, &consensus, 3, RecodeScheme::SR6);

        assert_eq!(rows[0], b"XXXXX".to_vec());
    }

    #[test]
    fn missing_ambiguous_and_invalid_positions_break_difference_runs() {
        let consensus = vec![Some(0); 9];
        let mut rows = vec![b"DDAD-DXBD".to_vec()];

        mask_rows_with_long_divergence(&mut rows, &consensus, 3, RecodeScheme::SR6);

        assert_eq!(rows[0], b"DDAD-DXBD".to_vec());
    }

    #[test]
    fn consensus_positions_without_valid_recoded_state_break_difference_runs() {
        let consensus = vec![Some(0), Some(0), None, Some(0), Some(0)];
        let mut rows = vec![b"DDDDA".to_vec()];

        mask_rows_with_long_divergence(&mut rows, &consensus, 3, RecodeScheme::SR6);

        assert_eq!(rows[0], b"DDDDA".to_vec());
    }

    #[test]
    fn mask_zero_returns_without_changing_rows() {
        let consensus = vec![Some(0); 5];
        let mut rows = vec![b"DDDDD".to_vec()];

        mask_rows_with_long_divergence(&mut rows, &consensus, 0, RecodeScheme::SR6);

        assert_eq!(rows[0], b"DDDDD".to_vec());
    }

    #[test]
    fn retained_middle_starts_at_first_retained_polymorphic_column() {
        let rows = &[b"ACCEF".as_slice(), b"ACCHF", b"ACCEF", b"ACCHF"];
        let result = filter_groups(sorted(rows), 1.0, 0, 0, 1).unwrap();

        assert_eq!(result.partitions, vec![(0, 0, 1)]);
        assert_eq!(result.concat[0], b"E".to_vec());
        assert_eq!(result.concat[1], b"H".to_vec());
        assert_eq!(result.concat[2], b"E".to_vec());
        assert_eq!(result.concat[3], b"H".to_vec());
    }

    #[test]
    fn retained_middle_trims_past_raw_only_polymorphism() {
        let rows = &[b"AAPCF".as_slice(), b"APDCF", b"AAPCF", b"APDCF"];
        let result = filter_groups(sorted(rows), 1.0, 0, 0, 1).unwrap();

        assert_eq!(result.partitions, vec![(0, 0, 1)]);
        assert_eq!(result.concat[0], b"P".to_vec());
        assert_eq!(result.concat[1], b"D".to_vec());
        assert_eq!(result.concat[2], b"P".to_vec());
        assert_eq!(result.concat[3], b"D".to_vec());
    }

    #[test]
    fn ambiguous_x_does_not_create_early_polymorphism() {
        let rows = &[b"ACCEF".as_slice(), b"ACCHF", b"AXCEF", b"ACXHF"];
        let result = filter_groups(sorted(rows), 0.75, 0, 0, 1).unwrap();

        assert_eq!(result.partitions, vec![(0, 0, 1)]);
        assert_eq!(result.concat[0], b"E".to_vec());
        assert_eq!(result.concat[1], b"H".to_vec());
        assert_eq!(result.concat[2], b"E".to_vec());
        assert_eq!(result.concat[3], b"H".to_vec());
    }

    #[test]
    fn discards_when_only_retained_middle_columns_are_conserved() {
        let rows = &[b"ACCEF".as_slice(), b"ACCDF", b"ACC-F", b"ACC-F"];
        let result = filter_groups(sorted(rows), 0.75, 0, 0, 1).unwrap();

        assert!(result.partitions.is_empty());
        assert!(result.concat.iter().all(Vec::is_empty));
    }

    #[test]
    fn right_anchor_tail_is_appended_after_polymorphic_middle() {
        let rows = &[b"ACEFT".as_slice(), b"ACHFT", b"ACEFT", b"ACHFT"];
        let result = filter_groups(sorted(rows), 1.0, 1, 0, 1).unwrap();

        assert_eq!(result.partitions, vec![(0, 1, 2)]);
        assert_eq!(result.concat[0], b"ET".to_vec());
        assert_eq!(result.concat[1], b"HT".to_vec());
        assert_eq!(result.concat[2], b"ET".to_vec());
        assert_eq!(result.concat[3], b"HT".to_vec());
    }

    #[test]
    fn polymorphic_right_anchor_tail_cannot_rescue_conserved_middle() {
        let rows = &[b"ACCD".as_slice(), b"ACCE", b"ACCD", b"ACCE"];
        let result = filter_groups(sorted(rows), 1.0, 1, 0, 1).unwrap();

        assert!(result.partitions.is_empty());
        assert!(result.concat.iter().all(Vec::is_empty));
    }

    #[test]
    fn discards_when_only_raw_polymorphic_middle_column_fails_missingness() {
        let rows = &[b"ACGT".as_slice(), b"ADGT", b"A-GT", b"A-GT"];
        let result = filter_groups(sorted(rows), 0.75, 0, 0, 1).unwrap();

        assert!(result.partitions.is_empty());
        assert!(result.concat.iter().all(Vec::is_empty));
    }

    #[test]
    fn retains_group_with_retained_polymorphic_middle_column() {
        let rows = &[b"ACGT".as_slice(), b"ADGT", b"ACGT", b"ADGT"];
        let result = filter_groups(sorted(rows), 0.75, 1, 0, 1).unwrap();

        assert_eq!(result.partitions, vec![(0, 1, 2)]);
        assert_eq!(result.concat[0], b"CT".to_vec());
        assert_eq!(result.concat[1], b"DT".to_vec());
        assert_eq!(result.concat[2], b"CT".to_vec());
        assert_eq!(result.concat[3], b"DT".to_vec());
    }

    #[test]
    fn constant_tail_cannot_rescue_without_retained_polymorphic_middle() {
        let rows = &[b"AC-T".as_slice(), b"AD-T", b"A--T", b"A--T"];
        let result = filter_groups(sorted(rows), 0.75, 1, 0, 1).unwrap();

        assert!(result.partitions.is_empty());
        assert!(result.concat.iter().all(Vec::is_empty));
    }

    #[test]
    fn middle_is_trimmed_on_both_sides_before_appending_constants() {
        let rows = &[b"AAAQFCPQRRR".as_slice(), b"AAAQDCDQRRR"];
        let result = filter_groups(sorted_with_k(rows, 3), 1.0, 3, 0, 1).unwrap();

        assert_eq!(result.partitions, vec![(0, 5, 6)]);
        assert_eq!(result.concat[0], b"FCPRRR".to_vec());
        assert_eq!(result.concat[1], b"DCDRRR".to_vec());
    }

    #[test]
    fn final_middle_column_before_constants_is_polymorphic() {
        let rows = &[b"AAAQFCPQRRR".as_slice(), b"AAAQDCDQRRR"];
        let result = filter_groups(sorted_with_k(rows, 3), 1.0, 3, 0, 1).unwrap();
        let tail_keep = 3;
        let middle_len = result.partitions[0].2 - tail_keep;
        let trimmed_rows: Vec<Vec<u8>> = result
            .concat
            .iter()
            .map(|row| row[..middle_len].to_vec())
            .collect();

        assert!(is_polymorphic_column(
            &trimmed_rows,
            middle_len - 1,
            RecodeScheme::SR6
        ));
        assert_eq!(result.concat[0][middle_len - 1], b'P');
        assert_eq!(result.concat[1][middle_len - 1], b'D');
    }

    #[test]
    fn missingness_filtering_happens_before_polymorphic_trimming() {
        let rows = &[b"AFCPR".as_slice(), b"ADCDR", b"A-CPR", b"A-CDR"];
        let result = filter_groups(sorted(rows), 0.75, 0, 0, 1).unwrap();

        assert_eq!(result.partitions, vec![(0, 0, 1)]);
        assert_eq!(result.concat[0], b"P".to_vec());
        assert_eq!(result.concat[1], b"D".to_vec());
        assert_eq!(result.concat[2], b"P".to_vec());
        assert_eq!(result.concat[3], b"D".to_vec());
    }

    #[test]
    fn raw_amino_acid_differences_in_same_recoded_state_do_not_count_as_polymorphism() {
        let rows = &[b"AACA".as_slice(), b"APCA"];
        let result = filter_groups(sorted(rows), 1.0, 0, 0, 1).unwrap();

        assert!(result.partitions.is_empty());
        assert!(result.concat.iter().all(Vec::is_empty));
    }

    #[test]
    fn gaps_and_x_are_ignored_when_detecting_polymorphism() {
        let rows = &[b"ACPR".as_slice(), b"AXPR", b"A-DR", b"ACDR"];
        let result = filter_groups(sorted(rows), 0.5, 0, 0, 1).unwrap();

        assert_eq!(result.partitions, vec![(0, 0, 1)]);
        assert_eq!(result.concat[0], b"P".to_vec());
        assert_eq!(result.concat[1], b"P".to_vec());
        assert_eq!(result.concat[2], b"D".to_vec());
        assert_eq!(result.concat[3], b"D".to_vec());
    }

    #[test]
    fn constant_zero_emits_only_trimmed_polymorphic_middle_region() {
        let rows = &[b"AAAQFCPQRRR".as_slice(), b"AAAQDCDQRRR"];
        let result = filter_groups(sorted_with_k(rows, 3), 1.0, 0, 0, 1).unwrap();

        assert_eq!(result.partitions, vec![(0, 2, 3)]);
        assert_eq!(result.concat[0], b"FCP".to_vec());
        assert_eq!(result.concat[1], b"DCD".to_vec());
    }

    #[test]
    fn one_polymorphic_middle_column_with_constant_three_has_length_four() {
        let rows = &[b"AAACRRR".as_slice(), b"AAADRRR"];
        let result = filter_groups(sorted_with_k(rows, 3), 1.0, 3, 0, 1).unwrap();

        assert_eq!(result.partitions, vec![(0, 3, 4)]);
        assert_eq!(result.concat[0], b"CRRR".to_vec());
        assert_eq!(result.concat[1], b"DRRR".to_vec());
    }

    #[test]
    fn divergent_row_masking_happens_before_polymorphic_middle_filtering() {
        let rows = &[b"ACDEFGK".as_slice(), b"ACDEFGK", b"APPPPPK"];
        let result = filter_groups(sorted(rows), 1.0, 0, 5, 1).unwrap();

        assert!(result.partitions.is_empty());
        assert!(result.concat.iter().all(Vec::is_empty));
    }

    #[test]
    fn mask_zero_disables_divergent_row_masking() {
        let rows = &[b"ACDEFGK".as_slice(), b"ACDEFGK", b"APPPPPK"];
        let result = filter_groups(sorted(rows), 1.0, 0, 0, 1).unwrap();

        assert_eq!(result.partitions, vec![(0, 4, 5)]);
        assert_eq!(result.concat[0], b"CDEFG".to_vec());
        assert_eq!(result.concat[1], b"CDEFG".to_vec());
        assert_eq!(result.concat[2], b"PPPPP".to_vec());
    }

    #[test]
    fn constant_tail_columns_failing_missingness_are_not_appended() {
        let rows = &[b"ACT".as_slice(), b"AD-", b"AC-", b"AD-"];
        let result = filter_groups(sorted(rows), 0.75, 1, 0, 1).unwrap();

        assert_eq!(result.partitions, vec![(0, 0, 1)]);
        assert_eq!(result.concat[0], b"C".to_vec());
        assert_eq!(result.concat[1], b"D".to_vec());
        assert_eq!(result.concat[2], b"C".to_vec());
        assert_eq!(result.concat[3], b"D".to_vec());
    }
}
