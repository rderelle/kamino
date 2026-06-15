//! Post-extraction raw group selection, sorting, and non-overlap filtering.
use hashbrown::HashSet;

use anyhow::Result;

use crate::group_extraction::{AnchorPair, BlockArena, BlockObs, RawExtractedGroups};
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
    selected: Vec<BlockObs>,
    len: usize,
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
    /// Full raw block length before consensus collapse and column filtering.
    pub(crate) raw_len: usize,
    /// Raw observations with the selected best-supported block length.
    pub(crate) selected: Vec<BlockObs>,
}

fn select_unique_best_supported_length(
    key: AnchorPair,
    obs: &[BlockObs],
    n_species: usize,
    min_needed: usize,
    k: usize,
) -> Result<LengthFilteredGroup, LengthSupportReject> {
    let mut by_len: Vec<LengthSupport> = Vec::with_capacity(4);
    for o in obs {
        let len = o.aa.len();
        if let Some(bucket) = by_len.iter_mut().find(|bucket| bucket.len == len) {
            bucket.species.insert(o.sid as usize);
        } else {
            let mut species = SpeciesMask::new(n_species);
            species.insert(o.sid as usize);
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
    if best.0 < 2 * k + 1 {
        return Err(LengthSupportReject::ShortBlock);
    }

    let selected = obs
        .iter()
        .filter(|o| o.aa.len() == best.0)
        .cloned()
        .collect();

    Ok(LengthFilteredGroup {
        key,
        selected,
        len: best.0,
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

fn path_has_used_kmer(
    observations: &[BlockObs],
    arena: &BlockArena,
    k: usize,
    recode_scheme: RecodeScheme,
    used_global: &HashSet<u64>,
) -> bool {
    // Stream every valid recoded k-mer represented by the raw selected block
    // observations and stop as soon as any k-mer has already been retained by a
    // stronger raw group. This preserves full-path overlap semantics without
    // allocating, sorting, or deduplicating a per-candidate k-mer vector.
    let mask = kmer_mask(k);

    for obs in observations {
        let aa = arena.get(obs.aa);

        if aa.len() < k {
            continue;
        }

        let mut val = 0u64;
        let mut have = 0usize;

        for &b in aa {
            let c = recode_byte(b, recode_scheme);

            if c == 255 {
                val = 0;
                have = 0;
                continue;
            }

            val = ((val << RECODE_BITS_PER_SYMBOL) | c as u64) & mask;

            if have < k {
                have += 1;
            }

            if have == k && used_global.contains(&val) {
                return true;
            }
        }
    }

    false
}

fn insert_path_kmers(
    observations: &[BlockObs],
    arena: &BlockArena,
    k: usize,
    recode_scheme: RecodeScheme,
    used_global: &mut HashSet<u64>,
) {
    // Insert valid rolling recoded k-mers directly into the global overlap set.
    // Duplicate k-mers are intentionally left to HashSet::insert so retained
    // candidates avoid the old temporary Vec sort/dedup hot path.
    let mask = kmer_mask(k);

    for obs in observations {
        let aa = arena.get(obs.aa);

        if aa.len() < k {
            continue;
        }

        let mut val = 0u64;
        let mut have = 0usize;

        for &b in aa {
            let c = recode_byte(b, recode_scheme);

            if c == 255 {
                val = 0;
                have = 0;
                continue;
            }

            val = ((val << RECODE_BITS_PER_SYMBOL) | c as u64) & mask;

            if have < k {
                have += 1;
            }

            if have == k {
                used_global.insert(val);
            }
        }
    }
}

fn retain_nonoverlapping_raw_candidates(
    mut raw_candidates: Vec<RawDirectCandidate>,
    arena: &BlockArena,
    overlap_k: usize,
    recode_scheme: RecodeScheme,
) -> Vec<RawDirectCandidate> {
    raw_candidates.sort_unstable_by(|a, b| {
        b.coverage
            .cmp(&a.coverage)
            .then_with(|| b.raw_len.cmp(&a.raw_len))
            .then_with(|| a.key.0.cmp(&b.key.0))
            .then_with(|| a.key.1.cmp(&b.key.1))
    });

    let mut used_global: HashSet<u64> = HashSet::new();
    let mut nonoverlap_raw_candidates = Vec::new();

    for cand in raw_candidates {
        if path_has_used_kmer(
            &cand.selected,
            arena,
            overlap_k,
            recode_scheme,
            &used_global,
        ) {
            continue;
        }

        insert_path_kmers(
            &cand.selected,
            arena,
            overlap_k,
            recode_scheme,
            &mut used_global,
        );

        nonoverlap_raw_candidates.push(cand);
    }

    nonoverlap_raw_candidates
}

pub(crate) fn consensus_rows_by_isolate(
    observations: &[BlockObs],
    arena: &BlockArena,
    n_species: usize,
    len: usize,
) -> Vec<Vec<u8>> {
    // Merge duplicate observations from the same species. Agreement keeps the amino
    // acid; disagreement becomes `X` for downstream consensus-row filters.
    let mut rows = vec![vec![b'-'; len]; n_species];
    for obs in observations {
        let sid = obs.sid as usize;
        for (cell, &aa) in rows[sid].iter_mut().zip(arena.get(obs.aa)).take(len) {
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
    /// Anchor length used to split left/right anchors from middle columns.
    pub(crate) k: usize,
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
        min_needed,
        n_species: n,
    } = raw;

    let mut raw_candidates = Vec::with_capacity(groups.len());
    // Stage 5: for each original anchor pair, keep one unambiguous raw block
    // length present in enough species. These raw candidates are sorted by
    // strength inside recoded k-mer non-overlap filtering.
    for (key, obs) in groups {
        if let Ok(group) = select_unique_best_supported_length(key, &obs, n, min_needed, k) {
            raw_candidates.push(RawDirectCandidate {
                key: group.key,
                coverage: group.coverage,
                raw_len: group.len,
                selected: group.selected,
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
        recode_scheme,
    );
    Ok(SortedGroups {
        species_names,
        arena,
        raw_candidates: nonoverlap_raw_candidates,
        protein_names,
        k,
        min_needed,
        n_species: n,
        recode_scheme,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::group_extraction::RawDirectGroups;

    fn obs(arena: &mut BlockArena, sid: u16, aa: &[u8]) -> BlockObs {
        BlockObs {
            sid,
            protein_id: sid as u32,
            aa: arena.push(aa).unwrap(),
        }
    }

    fn raw_extracted(rows: &[&[u8]], k: usize, key: AnchorPair) -> RawExtractedGroups {
        let mut groups = RawDirectGroups::new();
        let mut arena = BlockArena::new();
        groups.insert(
            key,
            rows.iter()
                .enumerate()
                .map(|(sid, aa)| obs(&mut arena, sid as u16, aa))
                .collect(),
        );

        RawExtractedGroups {
            arena,
            species_names: (0..rows.len())
                .map(|sid| format!("species_{sid}"))
                .collect(),
            groups,
            protein_names: vec!["protein".to_string(); rows.len()],
            k,
            min_needed: rows.len(),
            n_species: rows.len(),
        }
    }

    fn raw_candidate(
        arena: &mut BlockArena,
        key: AnchorPair,
        coverage: usize,
        raw_len: usize,
        rows: &[&[u8]],
    ) -> RawDirectCandidate {
        RawDirectCandidate {
            key,
            coverage,
            raw_len,
            selected: rows
                .iter()
                .enumerate()
                .map(|(sid, aa)| obs(arena, sid as u16, aa))
                .collect(),
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
            vec![
                obs(&mut arena, 0, shared_min_block),
                obs(&mut arena, 1, shared_min_block),
            ],
        );
        groups.insert(
            (30, 40),
            vec![
                obs(&mut arena, 0, shared_min_block),
                obs(&mut arena, 1, shared_min_block),
            ],
        );

        let sorted = sort_and_deduplicate_groups(
            RawExtractedGroups {
                arena,
                species_names: vec!["species_0".to_string(), "species_1".to_string()],
                groups,
                protein_names: vec!["protein_0".to_string(), "protein_1".to_string()],
                k: 9,
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
        let observations = vec![
            obs(&mut arena, 0, b"ACD"),
            obs(&mut arena, 0, b"ACD"),
            obs(&mut arena, 1, b"ACE"),
        ];

        let rows = consensus_rows_by_isolate(&observations, &arena, 2, 3);

        assert_eq!(rows, vec![b"ACD".to_vec(), b"ACE".to_vec()]);
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
            RecodeScheme::SR6,
        );
        let retained_keys: Vec<_> = retained.into_iter().map(|cand| cand.key).collect();

        assert_eq!(retained_keys, vec![(10, 11), (30, 31)]);
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
            retain_nonoverlapping_raw_candidates(raw_candidates, &arena, 5, RecodeScheme::SR6);
        let retained_keys: Vec<_> = retained.into_iter().map(|cand| cand.key).collect();

        assert_eq!(retained_keys, vec![(10, 11), (30, 31), (40, 41)]);
    }

    #[test]
    fn length_filter_rejects_ties_low_support_and_short_blocks() {
        let mut arena = BlockArena::new();
        let tie_obs = vec![
            obs(&mut arena, 0, b"ACDEF"),
            obs(&mut arena, 1, b"ACDEF"),
            obs(&mut arena, 2, b"ACDEFG"),
            obs(&mut arena, 3, b"ACDEFG"),
        ];
        assert!(matches!(
            select_unique_best_supported_length((1, 2), &tie_obs, 4, 2, 2),
            Err(LengthSupportReject::Tie)
        ));

        let low_support_obs = vec![obs(&mut arena, 0, b"ACDEF")];
        assert!(matches!(
            select_unique_best_supported_length((1, 2), &low_support_obs, 4, 2, 2),
            Err(LengthSupportReject::LowSupport)
        ));

        let short_obs = vec![obs(&mut arena, 0, b"ACD"), obs(&mut arena, 1, b"ACD")];
        assert!(matches!(
            select_unique_best_supported_length((1, 2), &short_obs, 2, 2, 2),
            Err(LengthSupportReject::ShortBlock)
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
        assert_eq!(sorted.raw_candidates[0].raw_len, rows[0].len());
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
            RecodeScheme::SR6,
        );
        let retained_keys: Vec<_> = retained.into_iter().map(|cand| cand.key).collect();

        assert_eq!(retained_keys, vec![(10, 20), (30, 40)]);
    }

    #[test]
    fn raw_direct_groups_already_group_exact_anchor_pairs() {
        let mut groups = RawDirectGroups::new();
        let mut arena = BlockArena::new();
        groups
            .entry((1, 2))
            .or_default()
            .push(obs(&mut arena, 0, b"ACDEF"));
        groups
            .entry((1, 2))
            .or_default()
            .push(obs(&mut arena, 1, b"ACDEF"));
        groups
            .entry((1, 3))
            .or_default()
            .push(obs(&mut arena, 2, b"ACDEF"));

        assert_eq!(groups.len(), 2);
        assert_eq!(groups.get(&(1, 2)).unwrap().len(), 2);
        assert_eq!(groups.get(&(1, 3)).unwrap().len(), 1);
    }
}
