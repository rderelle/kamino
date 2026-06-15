//! Raw variant group extraction pipeline.
//!
//! The pipeline recodes amino-acid FASTA records into a compact six-state alphabet,
//! finds `k`-mer anchors shared by enough species, extracts short variable
//! blocks between local bubble-boundary anchors, and merges raw observations across species.
use anyhow::{anyhow, bail, Result};
use hashbrown::HashMap;
use rayon::{prelude::*, ThreadPoolBuilder};
use seq_io::fasta::{Reader as FastaReader, Record};
use std::collections::BTreeMap;
use std::io::{self, Write};
use std::sync::{mpsc, Mutex};

use crate::io::{open_fasta, SpeciesInput};
use crate::proba_filter::{
    AtomicCountMinSketch, BloomFilter, CountMinSketch, BLOOM_BITS, CMS_DEPTH, CMS_WIDTH,
};
use crate::recode::{recode_byte, RECODE_BITS_PER_SYMBOL};
use crate::RecodeScheme;

const START_CLUSTER_SPAN: usize = 10;

struct ProteomeProgress {
    current: Mutex<usize>,
    total: usize,
}

impl ProteomeProgress {
    fn new(total: usize) -> Self {
        Self {
            current: Mutex::new(0),
            total,
        }
    }

    fn increment(&self) {
        let mut current = self.current.lock().expect("progress mutex poisoned");
        *current += 1;
        self.draw(*current);
    }

    fn draw(&self, current: usize) {
        const WIDTH: usize = 45;
        const BLUE: &str = "\x1b[34m";
        const RESET: &str = "\x1b[0m";

        let filled = current
            .saturating_mul(WIDTH)
            .checked_div(self.total)
            .unwrap_or(WIDTH);
        let empty = WIDTH.saturating_sub(filled);
        eprint!(
            "\r  {BLUE}[{}{}]{RESET} {}/{} proteome files",
            "=".repeat(filled),
            " ".repeat(empty),
            current,
            self.total
        );
        let _ = io::stderr().flush();
    }

    fn finish(&self) {
        eprintln!();
    }
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct BlockAa {
    start: usize,
    len: u16,
}

impl BlockAa {
    #[inline]
    pub(crate) fn len(&self) -> usize {
        self.len as usize
    }
}

#[derive(Clone, Default, Debug)]
pub(crate) struct BlockArena {
    data: Vec<u8>,
}

impl BlockArena {
    pub(crate) fn new() -> Self {
        Self { data: Vec::new() }
    }

    pub(crate) fn push(&mut self, slice: &[u8]) -> Result<BlockAa> {
        let start = self.data.len();
        if slice.len() > u16::MAX as usize {
            bail!("amino-acid block length {} exceeds u16::MAX", slice.len());
        }
        self.data.extend_from_slice(slice);
        Ok(BlockAa {
            start,
            len: slice.len() as u16,
        })
    }

    pub(crate) fn append(&mut self, other: &BlockArena) -> usize {
        let offset = self.data.len();
        self.data.extend_from_slice(&other.data);
        offset
    }

    #[inline]
    pub(crate) fn get(&self, block: BlockAa) -> &[u8] {
        let start = block.start;
        let end = start + block.len as usize;
        &self.data[start..end]
    }
}

/// Pair of left/right recoded anchors that brackets one candidate variable block.
pub(crate) type AnchorPair = (u64, u64);
#[derive(Clone, Copy, Debug)]
/// One observed amino-acid block between the same anchor pair in one species.
pub(crate) struct BlockObs {
    /// Zero-based species index; kept compact because it is stored for every hit.
    pub sid: u16,
    /// Compact index into the global protein-name table.
    pub protein_id: u32,
    /// Original amino-acid characters spanning left anchor, middle, and right anchor.
    pub aa: BlockAa,
}
/// All raw observations grouped by their bracketing anchor pair.
pub(crate) type RawDirectGroups = HashMap<AnchorPair, Vec<BlockObs>>;

#[derive(Clone, Copy, Debug)]
/// One rolling `k`-mer hit and the sequence interval it covers.
struct KmerHit {
    kmer: u64,
    start: usize,
    end: usize,
    occupancy: u32,
    shared: bool,
}

#[derive(Clone, Copy, Debug)]
struct StartCandidate {
    hit: KmerHit,
}

fn is_better_start(candidate: StartCandidate, best: StartCandidate) -> bool {
    candidate.hit.occupancy > best.hit.occupancy
        || (candidate.hit.occupancy == best.hit.occupancy && candidate.hit.start < best.hit.start)
}

fn select_best_starts_per_local_cluster(starts: &[StartCandidate]) -> Vec<KmerHit> {
    let mut selected = Vec::new();

    if starts.is_empty() {
        return selected;
    }

    let mut cluster_start_pos = starts[0].hit.start;
    let mut best = starts[0];

    for &candidate in &starts[1..] {
        if candidate.hit.start < cluster_start_pos + START_CLUSTER_SPAN {
            if is_better_start(candidate, best) {
                best = candidate;
            }
        } else {
            selected.push(best.hit);
            cluster_start_pos = candidate.hit.start;
            best = candidate;
        }
    }

    selected.push(best.hit);
    selected
}

#[allow(clippy::too_many_arguments)]
fn extract_bubble_pairs(
    aa: &[u8],
    start_anchors: &[KmerHit],
    end_anchors: &[KmerHit],
    length_middle: usize,
    sid: u16,
    protein_id: u32,
    groups: &mut RawDirectGroups,
    arena: &mut BlockArena,
) -> Result<()> {
    let mut first_end = 0usize;

    for &left in start_anchors {
        while first_end < end_anchors.len() && end_anchors[first_end].start <= left.end {
            first_end += 1;
        }

        let mut j = first_end;
        let mut best_right: Option<KmerHit> = None;

        while j < end_anchors.len() {
            let right = end_anchors[j];

            if right.start <= left.end {
                j += 1;
                continue;
            }

            let middle_len = right.start - left.end;

            if middle_len > length_middle {
                break;
            }

            let replace = match best_right {
                None => true,
                Some(best) => {
                    right.occupancy > best.occupancy
                        || (right.occupancy == best.occupancy && right.start < best.start)
                }
            };

            if replace {
                best_right = Some(right);
            }

            j += 1;
        }

        let Some(right) = best_right else {
            continue;
        };

        groups
            .entry((left.kmer, right.kmer))
            .or_default()
            .push(BlockObs {
                sid,
                protein_id,
                aa: arena.push(&aa[left.start..right.end])?,
            });
    }

    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn extract_from_valid_segment(
    aa: &[u8],
    recoded: &[u8],
    k: usize,
    k_mask: u64,
    cms: &CountMinSketch,
    min_needed: u32,
    length_middle: usize,
    sid: u16,
    protein_id: u32,
    groups: &mut RawDirectGroups,
    arena: &mut BlockArena,
) -> Result<()> {
    // This function works only on already-valid amino-acid stretches. Unknown or
    // ambiguous residues split proteins into separate valid segments upstream.
    if recoded.len() < k || k == 0 {
        return Ok(());
    }

    let mut start_candidates = Vec::new();
    let mut end_anchors = Vec::new();
    let mut previous_hit: Option<KmerHit> = None;

    // Roll a compact `k`-mer code over the segment and query the completed,
    // immutable CMS snapshot from the first pass. Streaming neighbour comparisons
    // are equivalent to the previous full-hit-vector boundary detection, but avoid
    // storing every hit in the segment.
    let mut roll = 0u64;
    let mut have = 0usize;
    for (pos, &a) in recoded.iter().enumerate() {
        debug_assert_ne!(a, 255);
        roll = ((roll << RECODE_BITS_PER_SYMBOL) | (a as u64)) & k_mask;
        if have < k {
            have += 1;
        }
        if have == k {
            let occupancy = cms.estimate(roll);
            let shared = occupancy >= min_needed;
            let start = pos + 1 - k;
            let current_hit = KmerHit {
                kmer: roll,
                start,
                end: start + k,
                occupancy,
                shared,
            };

            if let Some(previous) = previous_hit {
                if previous.shared && previous.occupancy > current_hit.occupancy {
                    start_candidates.push(StartCandidate { hit: previous });
                }
                if current_hit.shared && previous.occupancy < current_hit.occupancy {
                    end_anchors.push(current_hit);
                }
            }

            previous_hit = Some(current_hit);
        }
    }

    // Start candidates are shared k-mers preceding a lower-occupancy k-mer, ranked
    // by occupancy and then the left-most position within each local cluster. End
    // candidates are shared k-mers following a lower-occupancy k-mer, ranked by
    // occupancy and then the left-most position among valid right anchors. Starts
    // and ends are collected in hit order, so the extraction loop can prune by
    // position without sorting or testing every possible pair.
    let start_anchors = select_best_starts_per_local_cluster(&start_candidates);
    extract_bubble_pairs(
        aa,
        &start_anchors,
        &end_anchors,
        length_middle,
        sid,
        protein_id,
        groups,
        arena,
    )
}

struct SpeciesExtraction {
    /// Species index that lets parallel results be merged deterministically.
    sid: usize,
    /// Raw amino-acid bytes referenced by this species' block observations.
    arena: BlockArena,
    /// Raw block observations produced by this species.
    groups: RawDirectGroups,
    /// Protein descriptions indexed by per-observation protein IDs.
    protein_names: Vec<String>,
}

fn merge_species_extraction(
    extraction: SpeciesExtraction,
    groups: &mut RawDirectGroups,
    arena: &mut BlockArena,
    protein_names: &mut Vec<String>,
) -> Result<()> {
    let protein_offset = protein_names.len();
    if protein_offset + extraction.protein_names.len() > u32::MAX as usize {
        bail!("too many protein records across all species to index with u32");
    }
    protein_names.extend(extraction.protein_names);
    let arena_offset = arena.append(&extraction.arena);
    for (key, observations) in extraction.groups {
        let mut merged_observations = Vec::with_capacity(observations.len());
        for mut obs in observations {
            obs.protein_id += protein_offset as u32;
            obs.aa.start += arena_offset;
            merged_observations.push(obs);
        }
        groups.entry(key).or_default().extend(merged_observations);
    }
    Ok(())
}

fn count_species_kmers(
    input: &SpeciesInput,
    k: usize,
    k_mask: u64,
    recode_scheme: RecodeScheme,
    cmsa: &AtomicCountMinSketch,
) -> Result<()> {
    // Count each distinct `k`-mer at most once per species, which makes the
    // sketch estimate species occupancy rather than raw copy number.
    let rdr = open_fasta(&input.path)?;
    let mut r = FastaReader::new(rdr);
    let mut bloom = BloomFilter::new(BLOOM_BITS);
    while let Some(rec) = r.next() {
        let rec = rec?;
        let mut roll = 0u64;
        let mut have = 0usize;
        for &b in rec.seq() {
            if b.is_ascii_whitespace() {
                continue;
            }
            let a = recode_byte(b.to_ascii_uppercase(), recode_scheme);
            if a == 255 {
                roll = 0;
                have = 0;
                continue;
            }
            roll = ((roll << RECODE_BITS_PER_SYMBOL) | (a as u64)) & k_mask;
            if have < k {
                have += 1;
            }
            if have == k && !bloom.test_and_set(roll) {
                cmsa.increment(roll);
            }
        }
    }
    Ok(())
}

fn truncate_protein_name(raw: &str) -> String {
    let trimmed = raw.trim();
    trimmed
        .split_once(" [")
        .map(|(name, _)| name.trim_end())
        .unwrap_or(trimmed)
        .to_string()
}

#[allow(clippy::too_many_arguments)]
fn extract_species_blocks(
    sid: usize,
    input: &SpeciesInput,
    k: usize,
    k_mask: u64,
    cms: &CountMinSketch,
    min_needed: u32,
    length_middle: usize,
    recode_scheme: RecodeScheme,
) -> Result<SpeciesExtraction> {
    // Split each protein at ambiguous residues. Valid stretches are independently
    // scanned for adjacent shared-anchor runs.
    let rdr = open_fasta(&input.path)?;
    let mut r = FastaReader::new(rdr);
    let mut groups: RawDirectGroups = HashMap::new();
    let mut arena = BlockArena::new();
    let mut protein_names = Vec::new();
    let min_block_len = 2 * k + 1;
    let mut aa_segment = Vec::new();
    let mut recoded_segment = Vec::new();
    while let Some(rec) = r.next() {
        let rec = rec?;
        let protein_record_id = rec.id()?.trim();
        let protein_desc = rec.desc().transpose()?.map(str::trim).unwrap_or_default();
        let raw_protein_name = if protein_desc.is_empty() {
            protein_record_id.to_string()
        } else {
            format!("{protein_record_id} {protein_desc}")
        };
        let protein_name = truncate_protein_name(&raw_protein_name);
        let protein_id = protein_names.len();
        if protein_id > u32::MAX as usize {
            bail!(
                "too many protein records in species {} to index with u32",
                input.name
            );
        }
        protein_names.push(protein_name);
        let protein_id = protein_id as u32;
        aa_segment.clear();
        recoded_segment.clear();
        let seq_len = rec.seq().len();
        if aa_segment.capacity() < seq_len {
            aa_segment.reserve(seq_len - aa_segment.capacity());
        }
        if recoded_segment.capacity() < seq_len {
            recoded_segment.reserve(seq_len - recoded_segment.capacity());
        }
        for &b in rec.seq() {
            if b.is_ascii_whitespace() {
                continue;
            }
            let up = b.to_ascii_uppercase();
            let rv = recode_byte(up, recode_scheme);
            if rv == 255 {
                if recoded_segment.len() >= min_block_len {
                    extract_from_valid_segment(
                        &aa_segment,
                        &recoded_segment,
                        k,
                        k_mask,
                        cms,
                        min_needed,
                        length_middle,
                        sid as u16,
                        protein_id,
                        &mut groups,
                        &mut arena,
                    )?;
                }
                aa_segment.clear();
                recoded_segment.clear();
            } else {
                aa_segment.push(up);
                recoded_segment.push(rv);
            }
        }
        if recoded_segment.len() >= min_block_len {
            extract_from_valid_segment(
                &aa_segment,
                &recoded_segment,
                k,
                k_mask,
                cms,
                min_needed,
                length_middle,
                sid as u16,
                protein_id,
                &mut groups,
                &mut arena,
            )?;
        }
    }

    Ok(SpeciesExtraction {
        sid,
        arena,
        groups,
        protein_names,
    })
}

/// Handoff from raw group extraction to post-extraction sorting and deduplication.
pub(crate) struct RawExtractedGroups {
    /// Species names in deterministic row order.
    pub(crate) species_names: Vec<String>,
    /// Global amino-acid bytes referenced by raw observations.
    pub(crate) arena: BlockArena,
    /// Raw block observations grouped by original bracketing anchor pairs.
    pub(crate) groups: RawDirectGroups,
    /// Global protein-name table referenced by raw observations.
    pub(crate) protein_names: Vec<String>,
    /// Anchor length used to split left/right anchors from middle columns.
    pub(crate) k: usize,
    /// Minimum species support threshold used during extraction.
    pub(crate) min_needed: usize,
    /// Number of species represented by `species_names`.
    pub(crate) n_species: usize,
}
#[allow(clippy::too_many_arguments)]
pub(crate) fn extract_groups(
    inputs: &[SpeciesInput],
    k: usize,
    min_freq: f32,
    length_middle: usize,
    recode_scheme: RecodeScheme,
    num_threads: usize,
) -> Result<RawExtractedGroups> {
    // Stage 1: determine the occupancy threshold and build a deterministic local
    // Rayon pool instead of changing the process-global pool.
    let n = inputs.len();
    if n > u16::MAX as usize {
        bail!(
            "too many species ({n}) for compact direct block observations; maximum is {}",
            u16::MAX
        );
    }
    let min_needed = (min_freq * n.max(1) as f32).ceil() as u32;
    let k_mask = if k * RECODE_BITS_PER_SYMBOL as usize >= 64 {
        u64::MAX
    } else {
        (1u64 << (k * RECODE_BITS_PER_SYMBOL as usize)) - 1
    };
    let n_threads = num_threads.max(1);
    let pool = ThreadPoolBuilder::new().num_threads(n_threads).build()?;

    // Stage 2: parallel species-level counting into an atomic count-min sketch.
    eprintln!(" . pass 1: estimate k-mer occupancies");
    let count_progress = ProteomeProgress::new(n);
    count_progress.draw(0);
    let cmsa = AtomicCountMinSketch::new(CMS_WIDTH, CMS_DEPTH);
    let count_results = pool.install(|| {
        inputs
            .par_iter()
            .map(|input| {
                let result = count_species_kmers(input, k, k_mask, recode_scheme, &cmsa);
                count_progress.increment();
                result
            })
            .collect::<Vec<_>>()
    });
    count_progress.finish();
    for result in count_results {
        result?;
    }

    // Stage 3: snapshot immutable counts and extract candidate blocks per species.
    // Extraction intentionally uses this completed snapshot; counting and extraction
    // remain separate passes so anchors are evaluated against all species.
    let cms = cmsa.snapshot();

    // Stage 4: stream completed species extractions into one deterministic global
    // merge. Out-of-order species are buffered only until all lower sids are ready.
    eprintln!(" . pass 2: extract sequences");
    let extraction_progress = ProteomeProgress::new(n);
    extraction_progress.draw(0);
    let mut groups: RawDirectGroups = HashMap::new();
    let mut arena = BlockArena::new();
    let mut protein_names = Vec::new();
    std::thread::scope(|thread_scope| {
        let (tx, rx) = mpsc::channel();
        let extraction_thread = thread_scope.spawn(|| {
            pool.install(|| {
                inputs
                    .par_iter()
                    .enumerate()
                    .for_each_with(tx, |tx, (sid, input)| {
                        let result = extract_species_blocks(
                            sid,
                            input,
                            k,
                            k_mask,
                            &cms,
                            min_needed,
                            length_middle,
                            recode_scheme,
                        )
                        .map(|extraction| (extraction.sid, extraction));
                        extraction_progress.increment();
                        let _ = tx.send(result);
                    });
            });
        });

        let mut pending = BTreeMap::new();
        let mut next_sid_to_merge = 0usize;
        let mut first_error = None;
        for result in rx {
            match result {
                Ok((sid, extraction)) if first_error.is_none() => {
                    pending.insert(sid, extraction);
                    while let Some(extraction) = pending.remove(&next_sid_to_merge) {
                        if let Err(err) = merge_species_extraction(
                            extraction,
                            &mut groups,
                            &mut arena,
                            &mut protein_names,
                        ) {
                            first_error.get_or_insert(err);
                            break;
                        }
                        next_sid_to_merge += 1;
                    }
                }
                Ok((_sid, _extraction)) => {}
                Err(err) => {
                    first_error.get_or_insert(err);
                }
            }
        }
        extraction_thread
            .join()
            .expect("species extraction worker thread panicked");
        extraction_progress.finish();
        if let Some(err) = first_error {
            Err(err)
        } else if next_sid_to_merge != n {
            Err(anyhow!(
                "species extraction ended before all species were merged ({next_sid_to_merge}/{n})"
            ))
        } else {
            Ok(())
        }
    })?;
    Ok(RawExtractedGroups {
        species_names: inputs.iter().map(|x| x.name.clone()).collect(),
        arena,
        groups,
        protein_names,
        k,
        min_needed: min_needed as usize,
        n_species: n,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn truncate_protein_name_keeps_record_id_and_removes_organism_suffix() {
        assert_eq!(
            truncate_protein_name("WP_012345.1 hypothetical protein [Escherichia coli]"),
            "WP_012345.1 hypothetical protein"
        );
        assert_eq!(
            truncate_protein_name("WP_012345.1 hypothetical protein"),
            "WP_012345.1 hypothetical protein"
        );
    }

    #[test]
    fn block_arena_push_and_get_preserve_exact_bytes() {
        let mut arena = BlockArena::new();
        let first = arena.push(b"ACD").unwrap();
        let second = arena.push(b"EFGH").unwrap();

        assert_eq!(first.len(), 3);
        assert_eq!(second.len(), 4);
        assert_eq!(arena.get(first), b"ACD");
        assert_eq!(arena.get(second), b"EFGH");
    }

    fn extraction(sid: usize, protein: &str, key: AnchorPair, aa: &[u8]) -> SpeciesExtraction {
        let mut arena = BlockArena::new();
        let block = arena.push(aa).unwrap();
        let mut groups = RawDirectGroups::new();
        groups.insert(
            key,
            vec![BlockObs {
                sid: sid as u16,
                protein_id: 0,
                aa: block,
            }],
        );
        SpeciesExtraction {
            sid,
            arena,
            groups,
            protein_names: vec![protein.to_string()],
        }
    }

    #[test]
    fn pending_species_merge_is_deterministic_when_results_arrive_out_of_order() {
        let key = (7, 11);
        let mut groups = RawDirectGroups::new();
        let mut arena = BlockArena::new();
        let mut protein_names = Vec::new();
        let mut pending = BTreeMap::new();
        let mut next_sid_to_merge = 0usize;

        for extraction in [
            extraction(1, "protein_1", key, b"BBBB"),
            extraction(0, "protein_0", key, b"AAAA"),
            extraction(2, "protein_2", key, b"CCCC"),
        ] {
            pending.insert(extraction.sid, extraction);
            while let Some(extraction) = pending.remove(&next_sid_to_merge) {
                merge_species_extraction(extraction, &mut groups, &mut arena, &mut protein_names)
                    .unwrap();
                next_sid_to_merge += 1;
            }
        }

        assert_eq!(protein_names, vec!["protein_0", "protein_1", "protein_2"]);
        let observations = groups.get(&key).unwrap();
        assert_eq!(
            observations.iter().map(|obs| obs.sid).collect::<Vec<_>>(),
            vec![0, 1, 2]
        );
        assert_eq!(
            observations
                .iter()
                .map(|obs| obs.protein_id)
                .collect::<Vec<_>>(),
            vec![0, 1, 2]
        );
        assert_eq!(arena.get(observations[0].aa), b"AAAA");
        assert_eq!(arena.get(observations[1].aa), b"BBBB");
        assert_eq!(arena.get(observations[2].aa), b"CCCC");
    }
}
