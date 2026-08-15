//! Post-extraction raw group selection, sorting, and non-overlap filtering.
use hashbrown::HashSet;

use anyhow::Result;

use crate::group_extraction::{AnchorPair, BlockArena, RawExtractedGroups, RawGroup};
use crate::recode::{recode_byte, RECODE_BITS_PER_SYMBOL};
use crate::RecodeScheme;

const DEFAULT_DIRECT_OVERLAP_K: usize = 21;

fn direct_overlap_k(k: usize) -> usize {
    // The shortest raw group spans two k-mer anchors plus one middle
    // residue. Fixed recoded 21-mers work for k >= 11; for smaller k values the
    // overlap detector must use the whole minimum-length block so shortest valid
    // groups still contribute comparable k-mers.
    if k < 10 {
        2 * k + 1
    } else {
        DEFAULT_DIRECT_OVERLAP_K
    }
}

#[derive(Clone, Debug)]
enum SpeciesMask {
    Inline(u64),
    Inline512([u64; 8]),
    Dynamic(Vec<u64>),
}

impl SpeciesMask {
    fn new(n_species: usize) -> Self {
        if n_species <= 64 {
            Self::Inline(0)
        } else if n_species <= 512 {
            Self::Inline512([0; 8])
        } else {
            Self::Dynamic(vec![0; n_species.div_ceil(64)])
        }
    }

    #[inline]
    fn insert(&mut self, sid: usize) {
        match self {
            Self::Inline(mask) => *mask |= 1u64 << sid,
            Self::Inline512(words) => words[sid >> 6] |= 1u64 << (sid & 63),
            Self::Dynamic(words) => words[sid >> 6] |= 1u64 << (sid & 63),
        }
    }

    #[inline]
    fn count(&self) -> usize {
        match self {
            Self::Inline(mask) => mask.count_ones() as usize,
            Self::Inline512(words) => words.iter().map(|w| w.count_ones() as usize).sum(),
            Self::Dynamic(words) => words.iter().map(|w| w.count_ones() as usize).sum(),
        }
    }
}

struct LengthSupport {
    len: usize,
    species: SpeciesMask,
}

struct LengthFilteredGroup {
    key: AnchorPair,
    group: RawGroup,
    middle_len: usize,
    coverage: usize,
}

enum LengthSupportReject {
    Tie,
    LowSupport,
    ShortBlock,
}

pub(crate) struct RawDirectCandidate {
    /// Anchor pair used for deterministic ordering and early overlap checks.
    pub(crate) key: AnchorPair,
    /// Number of species supporting the selected block length.
    pub(crate) coverage: usize,
    /// Length of the variable middle between the two anchors.
    pub(crate) middle_len: usize,
    /// Deduplicated sequences and ordered occurrences for the selected middle length.
    pub(crate) group: RawGroup,
}

fn select_unique_best_supported_length(
    key: AnchorPair,
    mut group: RawGroup,
    n_species: usize,
    min_needed: usize,
    constant: usize,
) -> Result<LengthFilteredGroup, LengthSupportReject> {
    let mut by_len: Vec<LengthSupport> = Vec::with_capacity(4);
    for occurrence in &group.occurrences {
        let block = group.sequences[occurrence.sequence_id as usize];
        let len = block
            .len()
            .checked_sub(constant)
            .ok_or(LengthSupportReject::ShortBlock)?;
        if let Some(bucket) = by_len.iter_mut().find(|bucket| bucket.len == len) {
            bucket.species.insert(occurrence.sid as usize);
        } else {
            let mut species = SpeciesMask::new(n_species);
            species.insert(occurrence.sid as usize);
            by_len.push(LengthSupport { len, species });
        }
    }

    let mut best = (0usize, 0usize);
    let mut tie = false;
    for bucket in by_len {
        let c = bucket.species.count();
        if c > best.1 {
            best = (bucket.len, c);
            tie = false;
        } else if c == best.1 {
            tie = true;
        }
    }

    if tie {
        return Err(LengthSupportReject::Tie);
    }
    if best.1 < min_needed {
        return Err(LengthSupportReject::LowSupport);
    }
    if best.0 < 1 {
        return Err(LengthSupportReject::ShortBlock);
    }

    let mut remap = vec![u32::MAX; group.sequences.len()];
    let mut sequences = Vec::with_capacity(group.sequences.len());
    for (old_id, block) in group.sequences.into_iter().enumerate() {
        if block.len().checked_sub(constant) == Some(best.0) {
            let new_id = sequences.len() as u32;
            remap[old_id] = new_id;
            sequences.push(block);
        }
    }
    group.sequences = sequences;
    group.occurrences.retain_mut(|occurrence| {
        let new_id = remap[occurrence.sequence_id as usize];
        if new_id == u32::MAX {
            false
        } else {
            occurrence.sequence_id = new_id;
            true
        }
    });

    Ok(LengthFilteredGroup {
        key,
        group,
        middle_len: best.0,
        coverage: best.1,
    })
}

#[inline]
fn kmer_mask(k: usize) -> u64 {
    let bit_len = k * RECODE_BITS_PER_SYMBOL as usize;
    if bit_len >= 64 {
        u64::MAX
    } else {
        (1u64 << bit_len) - 1
    }
}

#[allow(clippy::too_many_arguments)]
fn path_has_used_kmer(
    sequences: &[crate::group_extraction::BlockAa],
    key: AnchorPair,
    arena: &BlockArena,
    k: usize,
    anchor_k: usize,
    middle_len: usize,
    recode_scheme: RecodeScheme,
    used_global: &HashSet<u64>,
) -> bool {
    // Stream every valid recoded k-mer represented by the raw selected block
    // observations and stop as soon as any k-mer has already been retained by a
    // stronger raw group. This preserves full-path overlap semantics without
    // allocating, sorting, or deduplicating a per-candidate k-mer vector.
    let mask = kmer_mask(k);

    for &block in sequences {
        let mut val = 0u64;
        let mut have = 0usize;
        let mut feed = |c: u8| {
            if c == 255 {
                val = 0;
                have = 0;
                return false;
            }
            val = ((val << RECODE_BITS_PER_SYMBOL) | c as u64) & mask;
            if have < k {
                have += 1;
            }
            have == k && used_global.contains(&val)
        };
        for i in 0..anchor_k {
            if feed(((key.0 >> ((anchor_k - 1 - i) * RECODE_BITS_PER_SYMBOL as usize)) & 7) as u8) {
                return true;
            }
        }
        for &b in &arena.get(block)[..middle_len] {
            if feed(recode_byte(b, recode_scheme)) {
                return true;
            }
        }
        for i in 0..anchor_k {
            if feed(((key.1 >> ((anchor_k - 1 - i) * RECODE_BITS_PER_SYMBOL as usize)) & 7) as u8) {
                return true;
            }
        }
    }

    false
}

#[allow(clippy::too_many_arguments)]
fn insert_path_kmers(
    sequences: &[crate::group_extraction::BlockAa],
    key: AnchorPair,
    arena: &BlockArena,
    k: usize,
    anchor_k: usize,
    middle_len: usize,
    recode_scheme: RecodeScheme,
    used_global: &mut HashSet<u64>,
) {
    // Insert valid rolling recoded k-mers directly into the global overlap set.
    // Duplicate k-mers are intentionally left to HashSet::insert so retained
    // candidates avoid the old temporary Vec sort/dedup hot path.
    let mask = kmer_mask(k);

    for &block in sequences {
        let mut val = 0u64;
        let mut have = 0usize;
        let mut feed = |c: u8| {
            if c == 255 {
                val = 0;
                have = 0;
                return;
            }
            val = ((val << RECODE_BITS_PER_SYMBOL) | c as u64) & mask;
            if have < k {
                have += 1;
            }
            if have == k {
                used_global.insert(val);
            }
        };
        for i in 0..anchor_k {
            feed(((key.0 >> ((anchor_k - 1 - i) * RECODE_BITS_PER_SYMBOL as usize)) & 7) as u8);
        }
        for &b in &arena.get(block)[..middle_len] {
            feed(recode_byte(b, recode_scheme));
        }
        for i in 0..anchor_k {
            feed(((key.1 >> ((anchor_k - 1 - i) * RECODE_BITS_PER_SYMBOL as usize)) & 7) as u8);
        }
    }
}

fn retain_nonoverlapping_raw_candidates(
    mut raw_candidates: Vec<RawDirectCandidate>,
    arena: &BlockArena,
    overlap_k: usize,
    anchor_k: usize,
    recode_scheme: RecodeScheme,
) -> Vec<RawDirectCandidate> {
    raw_candidates.sort_unstable_by(|a, b| {
        b.coverage
            .cmp(&a.coverage)
            .then_with(|| b.middle_len.cmp(&a.middle_len))
            .then_with(|| a.key.0.cmp(&b.key.0))
            .then_with(|| a.key.1.cmp(&b.key.1))
    });

    let mut used_global: HashSet<u64> = HashSet::new();
    let mut nonoverlap_raw_candidates = Vec::new();

    for cand in raw_candidates {
        if path_has_used_kmer(
            &cand.group.sequences,
            cand.key,
            arena,
            overlap_k,
            anchor_k,
            cand.middle_len,
            recode_scheme,
            &used_global,
        ) {
            continue;
        }

        insert_path_kmers(
            &cand.group.sequences,
            cand.key,
            arena,
            overlap_k,
            anchor_k,
            cand.middle_len,
            recode_scheme,
            &mut used_global,
        );

        nonoverlap_raw_candidates.push(cand);
    }

    nonoverlap_raw_candidates
}

pub(crate) fn consensus_rows_by_isolate(
    group: &RawGroup,
    arena: &BlockArena,
    n_species: usize,
    len: usize,
) -> Vec<Vec<u8>> {
    // Merge duplicate observations from the same species. Agreement keeps the amino
    // acid; disagreement becomes `X` for downstream consensus-row filters.
    let mut rows = vec![vec![b'-'; len]; n_species];
    for occurrence in &group.occurrences {
        let sid = occurrence.sid as usize;
        let block = group.sequences[occurrence.sequence_id as usize];
        for (cell, &aa) in rows[sid].iter_mut().zip(arena.get(block)).take(len) {
            match *cell {
                b'-' => *cell = aa,
                old if old == aa => {}
                _ => *cell = b'X',
            }
        }
    }
    rows
}

pub(crate) struct SortedGroups {
    /// Species names in deterministic row order.
    pub(crate) species_names: Vec<String>,
    /// Global amino-acid bytes referenced by raw observations.
    pub(crate) arena: BlockArena,
    /// Deduplicated raw candidates retained after recoded k-mer non-overlap filtering.
    pub(crate) raw_candidates: Vec<RawDirectCandidate>,
    /// Global protein-name table referenced by raw observations.
    pub(crate) protein_names: Vec<String>,
    /// Length of each complete recoded anchor in the logical full path.
    pub(crate) k: usize,
    /// Stored right-anchor prefix length.
    pub(crate) constant: usize,
    /// Minimum species support threshold used during extraction.
    pub(crate) min_needed: usize,
    /// Number of species represented by `species_names`.
    pub(crate) n_species: usize,
    /// Recoding scheme used for anchor encoding and recoded-space polymorphism checks.
    pub(crate) recode_scheme: RecodeScheme,
}

pub(crate) fn sort_and_deduplicate_groups(
    raw: RawExtractedGroups,
    recode_scheme: RecodeScheme,
) -> Result<SortedGroups> {
    let RawExtractedGroups {
        species_names,
        arena,
        groups,
        protein_names,
        k,
        constant,
        min_needed,
        n_species: n,
    } = raw;

    let mut raw_candidates = Vec::with_capacity(groups.len());
    // Stage 5: for each original anchor pair, keep one unambiguous raw block
    // length present in enough species. These raw candidates are sorted by
    // strength inside recoded k-mer non-overlap filtering.
    for (key, group) in groups {
        if let Ok(group) = select_unique_best_supported_length(key, group, n, min_needed, constant)
        {
            raw_candidates.push(RawDirectCandidate {
                key: group.key,
                coverage: group.coverage,
                middle_len: group.middle_len,
                group: group.group,
            });
        }
    }

    // Stage 6: greedily keep the strongest non-overlapping raw groups so one
    // biological signal cannot contribute the same transient recoded k-mer
    // twice. The detector uses recoded 21-mers except when k < 11, where the
    // minimum valid raw block is shorter than 21 residues and therefore defines
    // the comparison k-mer length. Raw groups without valid overlap k-mers pass
    // because they cannot be compared by the overlap set.
    let nonoverlap_raw_candidates = retain_nonoverlapping_raw_candidates(
        raw_candidates,
        &arena,
        direct_overlap_k(k),
        k,
        recode_scheme,
    );
    Ok(SortedGroups {
        species_names,
        arena,
        raw_candidates: nonoverlap_raw_candidates,
        protein_names,
        k,
        constant,
        min_needed,
        n_species: n,
        recode_scheme,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::group_extraction::{BlockOccurrence, RawDirectGroups};

    fn group(arena: &mut BlockArena, rows: &[(u16, &[u8])]) -> RawGroup {
        let mut group = RawGroup::default();
        for &(sid, aa) in rows {
            let sequence_id = group
                .sequences
                .iter()
                .position(|&block| arena.get(block) == aa)
                .unwrap_or_else(|| {
                    group.sequences.push(arena.push(aa).unwrap());
                    group.sequences.len() - 1
                }) as u32;
            group.occurrences.push(BlockOccurrence {
                sid,
                protein_id: sid as u32,
                sequence_id,
            });
        }
        group
    }

    fn raw_extracted(rows: &[&[u8]], k: usize, key: AnchorPair) -> RawExtractedGroups {
        let mut groups = RawDirectGroups::new();
        let mut arena = BlockArena::new();
        groups.insert(
            key,
            group(
                &mut arena,
                &rows
                    .iter()
                    .enumerate()
                    .map(|(sid, &aa)| (sid as u16, aa))
                    .collect::<Vec<_>>(),
            ),
        );

        RawExtractedGroups {
            arena,
            species_names: (0..rows.len())
                .map(|sid| format!("species_{sid}"))
                .collect(),
            groups,
            protein_names: vec!["protein".to_string(); rows.len()],
            k,
            constant: 0,
            min_needed: rows.len(),
            n_species: rows.len(),
        }
    }

    fn raw_candidate(
        arena: &mut BlockArena,
        key: AnchorPair,
        coverage: usize,
        middle_len: usize,
        rows: &[&[u8]],
    ) -> RawDirectCandidate {
        RawDirectCandidate {
            key,
            coverage,
            middle_len,
            group: group(
                arena,
                &rows
                    .iter()
                    .enumerate()
                    .map(|(sid, &aa)| (sid as u16, aa))
                    .collect::<Vec<_>>(),
            ),
        }
    }

    #[test]
    fn direct_overlap_k_uses_minimum_block_length_below_k_11() {
        assert_eq!(direct_overlap_k(1), 3);
        assert_eq!(direct_overlap_k(9), 19);
        assert_eq!(direct_overlap_k(10), DEFAULT_DIRECT_OVERLAP_K);
        assert_eq!(direct_overlap_k(20), DEFAULT_DIRECT_OVERLAP_K);
    }

    #[test]
    fn species_mask_uses_inline_u64_for_up_to_64_species() {
        let mut mask = SpeciesMask::new(64);
        assert!(matches!(mask, SpeciesMask::Inline(_)));

        mask.insert(0);
        mask.insert(63);
        mask.insert(63);

        assert_eq!(mask.count(), 2);
    }

    #[test]
    fn species_mask_uses_inline512_for_65_to_512_species() {
        let mut mask = SpeciesMask::new(65);
        assert!(matches!(mask, SpeciesMask::Inline512(_)));
        mask.insert(0);
        mask.insert(64);
        assert_eq!(mask.count(), 2);

        let mut mask = SpeciesMask::new(512);
        assert!(matches!(mask, SpeciesMask::Inline512(_)));
        mask.insert(0);
        mask.insert(63);
        mask.insert(64);
        mask.insert(511);
        assert_eq!(mask.count(), 4);
    }

    #[test]
    fn species_mask_uses_dynamic_for_more_than_512_species() {
        let mut mask = SpeciesMask::new(513);
        assert!(matches!(mask, SpeciesMask::Dynamic(_)));

        mask.insert(0);
        mask.insert(512);
        mask.insert(512);

        assert_eq!(mask.count(), 2);
    }

    #[test]
    fn sort_nonoverlap_uses_shorter_overlap_k_when_k_less_than_11() {
        let mut groups = RawDirectGroups::new();
        let mut arena = BlockArena::new();
        let shared_min_block = b"ACDEFGHIKLMNPQRSTVW";
        groups.insert(
            (10, 20),
            group(&mut arena, &[(0, shared_min_block), (1, shared_min_block)]),
        );
        groups.insert(
            (30, 40),
            group(&mut arena, &[(0, shared_min_block), (1, shared_min_block)]),
        );

        let sorted = sort_and_deduplicate_groups(
            RawExtractedGroups {
                arena,
                species_names: vec!["species_0".to_string(), "species_1".to_string()],
                groups,
                protein_names: vec!["protein_0".to_string(), "protein_1".to_string()],
                k: 9,
                constant: 0,
                min_needed: 2,
                n_species: 2,
            },
            RecodeScheme::SR6,
        )
        .unwrap();

        let retained_keys: Vec<_> = sorted.raw_candidates.iter().map(|cand| cand.key).collect();
        assert_eq!(retained_keys, vec![(10, 20)]);
    }

    #[test]
    fn consensus_rows_by_isolate_uses_arena_backed_observations() {
        let mut arena = BlockArena::new();
        let observations = group(
            &mut arena,
            &[
                (0, b"ACD"),
                (1, b"ACD"),
                (1, b"ACD"),
                (2, b"ACD"),
                (2, b"AQD"),
            ],
        );

        let rows = consensus_rows_by_isolate(&observations, &arena, 4, 3);

        assert_eq!(observations.sequences.len(), 2);
        assert_eq!(observations.occurrences.len(), 5);
        assert_eq!(
            rows,
            vec![
                b"ACD".to_vec(),
                b"ACD".to_vec(),
                b"AXD".to_vec(),
                b"---".to_vec()
            ]
        );
    }

    #[test]
    fn nonoverlap_filter_drops_weaker_shared_kmer_and_retains_empty_kmer_candidates() {
        let mut arena = BlockArena::new();
        let shared_21mer = b"ACDEFGHIKLMNPQRSTVWYA";
        let raw_candidates = vec![
            raw_candidate(&mut arena, (20, 21), 1, shared_21mer.len(), &[shared_21mer]),
            raw_candidate(&mut arena, (30, 31), 1, 21, &[b"XXXXXXXXXXXXXXXXXXXXX"]),
            raw_candidate(&mut arena, (10, 11), 2, shared_21mer.len(), &[shared_21mer]),
        ];

        let retained = retain_nonoverlapping_raw_candidates(
            raw_candidates,
            &arena,
            DEFAULT_DIRECT_OVERLAP_K,
            0,
            RecodeScheme::SR6,
        );
        let retained_keys: Vec<_> = retained.into_iter().map(|cand| cand.key).collect();

        assert_eq!(retained_keys, vec![(10, 11), (30, 31)]);
    }

    #[test]
    fn nonoverlap_scans_one_unique_path_for_ten_identical_occurrences() {
        let mut arena = BlockArena::new();
        let sequence = b"ACDEFGHIKLMNPQRSTVWYA".as_slice();
        let rows = vec![sequence; 10];
        let candidate = raw_candidate(&mut arena, (1, 2), 10, sequence.len(), &rows);
        assert_eq!(candidate.group.sequences.len(), 1);
        assert_eq!(candidate.group.occurrences.len(), 10);
        let mut used = HashSet::new();
        insert_path_kmers(
            &candidate.group.sequences,
            candidate.key,
            &arena,
            DEFAULT_DIRECT_OVERLAP_K,
            0,
            candidate.middle_len,
            RecodeScheme::SR6,
            &mut used,
        );
        assert_eq!(used.len(), 1);
        assert!(path_has_used_kmer(
            &candidate.group.sequences,
            candidate.key,
            &arena,
            DEFAULT_DIRECT_OVERLAP_K,
            0,
            candidate.middle_len,
            RecodeScheme::SR6,
            &used,
        ));
    }

    #[test]
    fn streaming_nonoverlap_preserves_full_path_overlap_semantics() {
        let mut arena = BlockArena::new();
        let shared = b"ACDEF";
        let raw_candidates = vec![
            raw_candidate(&mut arena, (20, 21), 3, 9, &[b"ZZACDEFZZ"]),
            raw_candidate(&mut arena, (40, 41), 1, 5, &[b"GGGGG"]),
            raw_candidate(&mut arena, (30, 31), 1, 5, &[b"XXXXX"]),
            raw_candidate(&mut arena, (10, 11), 4, shared.len(), &[shared]),
        ];

        let retained =
            retain_nonoverlapping_raw_candidates(raw_candidates, &arena, 5, 0, RecodeScheme::SR6);
        let retained_keys: Vec<_> = retained.into_iter().map(|cand| cand.key).collect();

        assert_eq!(retained_keys, vec![(10, 11), (30, 31), (40, 41)]);
    }

    #[test]
    fn compact_path_reconstruction_matches_historical_full_path_kmers() {
        fn encode(aa: &[u8]) -> u64 {
            aa.iter().fold(0, |v, &b| {
                (v << RECODE_BITS_PER_SYMBOL) | recode_byte(b, RecodeScheme::SR6) as u64
            })
        }
        let left = b"ACD";
        let middle = b"EF";
        let right = b"GHI";
        let key = (encode(left), encode(right));
        let mut arena = BlockArena::new();
        // The stored I is right-anchor context and must not be streamed twice.
        let observations = group(&mut arena, &[(0, b"EFI")]);
        let mut reconstructed = HashSet::new();
        insert_path_kmers(
            &observations.sequences,
            key,
            &arena,
            3,
            3,
            middle.len(),
            RecodeScheme::SR6,
            &mut reconstructed,
        );

        let full = [left.as_slice(), middle.as_slice(), right.as_slice()].concat();
        let expected: HashSet<u64> = full.windows(3).map(encode).collect();
        assert_eq!(reconstructed, expected);
        // These include an anchor-only window and both anchor/middle boundaries.
        assert!(reconstructed.contains(&encode(left)));
        assert!(reconstructed.contains(&encode(b"DEF")));
        assert!(reconstructed.contains(&encode(b"FGH")));
        assert!(reconstructed.contains(&encode(right)));
    }

    #[test]
    fn length_filter_rejects_ties_low_support_and_short_blocks() {
        let mut arena = BlockArena::new();
        let tie_obs = group(
            &mut arena,
            &[(0, b"ACDEF"), (1, b"ACDEF"), (2, b"ACDEFG"), (3, b"ACDEFG")],
        );
        assert!(matches!(
            select_unique_best_supported_length((1, 2), tie_obs, 4, 2, 2),
            Err(LengthSupportReject::Tie)
        ));

        let low_support_obs = group(&mut arena, &[(0, b"ACDEF")]);
        assert!(matches!(
            select_unique_best_supported_length((1, 2), low_support_obs, 4, 2, 2),
            Err(LengthSupportReject::LowSupport)
        ));

        let short_obs = group(&mut arena, &[(0, b""), (1, b"")]);
        assert!(matches!(
            select_unique_best_supported_length((1, 2), short_obs, 2, 2, 0),
            Err(LengthSupportReject::ShortBlock)
        ));
    }

    #[test]
    fn length_support_counts_distinct_species_not_occurrences_or_sequences() {
        let mut arena = BlockArena::new();
        let observations = group(
            &mut arena,
            &[
                (0, b"AAAA"),
                (0, b"BBBB"),
                (0, b"AAAA"),
                (1, b"CCCC"),
                (2, b"DDDDD"),
                (3, b"EEEEE"),
            ],
        );
        assert!(matches!(
            select_unique_best_supported_length((1, 2), observations, 4, 2, 0),
            Err(LengthSupportReject::Tie)
        ));
    }

    #[test]
    fn sort_retains_original_raw_anchor_key_without_refinement() {
        let rows = &[b"*AC**FGHIKL*".as_slice(), b"*AC**DGHIKL*"];
        let sorted =
            sort_and_deduplicate_groups(raw_extracted(rows, 2, (31, 37)), RecodeScheme::SR6)
                .unwrap();

        assert_eq!(sorted.raw_candidates.len(), 1);
        assert_eq!(sorted.raw_candidates[0].key, (31, 37));
        assert_eq!(sorted.raw_candidates[0].middle_len, rows[0].len());
    }

    #[test]
    fn sort_nonoverlap_keeps_strongest_overlapping_group_deterministically() {
        let mut arena = BlockArena::new();
        let shared_21mer = b"ACDEFGHIKLMNPQRSTVWYA";
        let unique_21mer = b"YYYYYYYYYYYYYYYYYYYYY";
        let raw_candidates = vec![
            raw_candidate(&mut arena, (50, 60), 2, shared_21mer.len(), &[shared_21mer]),
            raw_candidate(&mut arena, (10, 20), 2, shared_21mer.len(), &[shared_21mer]),
            raw_candidate(&mut arena, (30, 40), 1, unique_21mer.len(), &[unique_21mer]),
        ];

        let retained = retain_nonoverlapping_raw_candidates(
            raw_candidates,
            &arena,
            DEFAULT_DIRECT_OVERLAP_K,
            0,
            RecodeScheme::SR6,
        );
        let retained_keys: Vec<_> = retained.into_iter().map(|cand| cand.key).collect();

        assert_eq!(retained_keys, vec![(10, 20), (30, 40)]);
    }

    #[test]
    fn raw_direct_groups_already_group_exact_anchor_pairs() {
        let mut groups = RawDirectGroups::new();
        let mut arena = BlockArena::new();
        groups.insert((1, 2), group(&mut arena, &[(0, b"ACDEF"), (1, b"ACDEF")]));
        groups.insert((1, 3), group(&mut arena, &[(2, b"ACDEF")]));

        assert_eq!(groups.len(), 2);
        assert_eq!(groups.get(&(1, 2)).unwrap().occurrences.len(), 2);
        assert_eq!(groups.get(&(1, 3)).unwrap().occurrences.len(), 1);
    }
}
