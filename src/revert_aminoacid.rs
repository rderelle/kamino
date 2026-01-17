use anyhow::{bail, Result};
use bitvec::prelude::*;
use hashbrown::HashMap;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use seq_io::fasta::{Reader as FastaReader, Record};
use std::path::PathBuf;

use crate::graph::Graph;
use crate::io::{open_fasta, SpeciesInput};
use crate::recode::{recode_byte, RecodeScheme, RECODE_BITS_PER_SYMBOL};
use crate::traverse::{PathRec, VariantGroups};
use std::sync::{Arc, Mutex};

pub(crate) type SpeciesKmerConsensus = (
    usize,
    Vec<HashMap<u64, usize>>,
    Vec<Vec<Option<Vec<u8>>>>,
    Vec<Vec<KmerNameStats>>,
    Vec<String>,
);

#[derive(Clone, Debug, Default)]
pub(crate) struct KmerNameStats {
    pub total_sets: u32,
    pub word_counts: Vec<(u32, u32)>,
}

impl KmerNameStats {
    fn add_words(&mut self, word_ids: &[u32]) {
        if word_ids.is_empty() {
            return;
        }
        self.total_sets = self.total_sets.saturating_add(1);
        for &word_id in word_ids {
            if let Some(entry) = self
                .word_counts
                .iter_mut()
                .find(|(existing_id, _)| *existing_id == word_id)
            {
                entry.1 = entry.1.saturating_add(1);
            } else {
                self.word_counts.push((word_id, 1));
            }
        }
    }
}

struct WordInterner {
    map: HashMap<String, u32>,
    words: Vec<String>,
}

impl WordInterner {
    fn new() -> Self {
        Self {
            map: HashMap::new(),
            words: Vec::new(),
        }
    }

    fn intern(&mut self, word: &str) -> u32 {
        if let Some(id) = self.map.get(word) {
            return *id;
        }
        let id = self.words.len() as u32;
        self.words.push(word.to_string());
        self.map.insert(self.words[id as usize].clone(), id);
        id
    }
}

enum PathData<'a> {
    Recoded {
        codes: Vec<u8>,
        species: &'a BitVec<u32, Lsb0>,
    },
}

/// Decode packed (k-1)-mer node into AA codes (MSB..LSB).
fn decode_k1_codes(node: u64, k: usize, sym_bits: u32) -> Vec<u8> {
    let k1 = k - 1;
    let mut out = vec![0u8; k1];
    for (i, slot) in out.iter_mut().enumerate().take(k1) {
        let shift = sym_bits * (k1 as u32 - 1 - i as u32);
        let code = ((node >> shift) & ((1u64 << sym_bits) - 1)) as u8;
        *slot = code;
    }
    out
}

/// Appended edge symbol is the low `sym_bits` of the successor node.
#[inline]
fn edge_sym_from_suffix(v: u64, sym_bits: u32) -> u32 {
    (v & ((1u64 << sym_bits) - 1)) as u32
}

/// Render path → AA/recoded codes: first node’s (k-1) + one symbol per edge.
fn render_path_codes(path: &[u64], k: usize, sym_bits: u32) -> Vec<u8> {
    debug_assert!(!path.is_empty());
    let mut codes: Vec<u8> = Vec::with_capacity((k - 1) + path.len() - 1);
    let k1_codes = decode_k1_codes(path[0], k, sym_bits);
    codes.extend_from_slice(&k1_codes);
    for w in path.windows(2) {
        let v = w[1];
        let sym = edge_sym_from_suffix(v, sym_bits) as u8;
        codes.push(sym);
    }
    codes
}

/// Return the mask covering the low `k * sym_bits` bits, saturating at 128 bits.
#[inline]
fn kmer_mask(sym_bits: u32, k: usize) -> u64 {
    let kb = sym_bits.saturating_mul(k as u32);
    if kb >= 64 {
        u64::MAX
    } else {
        (1u64 << kb) - 1
    }
}

/// Given a recoded symbol sequence `codes`, extract all length-k k-mers,
/// encoded as packed integers with `sym_bits` bits per symbol.
/// Unknown code (255) breaks the window (k-mers spanning it are skipped).
fn encoded_kmers_from_codes(codes: &[u8], k: usize, sym_bits: u32) -> Vec<u64> {
    let mut out = Vec::new();
    if k == 0 || codes.len() < k {
        return out;
    }

    out.reserve(codes.len() - k + 1);

    let mask = kmer_mask(sym_bits, k);

    let mut val: u64 = 0;
    let mut have: usize = 0;

    for &c in codes {
        if c == 255 {
            // Unknown symbol: reset the window.
            val = 0;
            have = 0;
            continue;
        }

        val = (val << sym_bits) | c as u64;

        if have + 1 < k {
            have += 1;
            continue;
        }

        if have + 1 > k {
            // Maintain only the last k symbols.
            val &= mask;
        } else {
            have = k;
        }

        out.push(val);
    }

    out
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn build_species_kmer_consensus(
    inputs: &[SpeciesInput],
    g: &Graph,
    groups: &VariantGroups,
    k: usize,
    head_max: usize,
    num_threads: usize,
) -> Result<SpeciesKmerConsensus> {
    let _ = head_max;
    let species = &g.species_names;
    let n = g.n_species;
    let sym_bits = RECODE_BITS_PER_SYMBOL;
    let recode_scheme: RecodeScheme = g.recode_scheme;

    // Set scan_k (=21 or less if short k-mer size)
    let min_vg_len: usize = 2 * (k - 1) + 1;
    let scan_k: usize = min_vg_len.min(21);

    // --------------------------------------------------------------------
    // PHASE 1: from variant groups, collect per-species k-mers
    // --------------------------------------------------------------------

    use hashbrown::hash_map::Entry;

    // species_kmer_maps[sid]: recoded k-mer (u64) -> index in species_kmer_consensus[sid]
    let mut species_kmer_maps: Vec<HashMap<u64, usize>> = vec![HashMap::new(); n];
    // species_kmer_consensus[sid][idx]: Option<Vec<u8>> = consensus 20AA k-mer
    let mut species_kmer_consensus: Vec<Vec<Option<Vec<u8>>>> = vec![Vec::new(); n];
    // species_kmer_name_stats[sid][idx]: KmerNameStats for protein name words
    let mut species_kmer_name_stats: Vec<Vec<KmerNameStats>> = vec![Vec::new(); n];

    for paths in groups.values() {
        for p in paths {
            // Recoded symbol sequence for this path (Dayhoff6 codes).
            let path_codes = render_path_codes(&p.nodes, k, sym_bits);
            let kmers = encoded_kmers_from_codes(&path_codes, scan_k, sym_bits);

            // For every species present on this path, register these k-mers.
            for (sid, ((map, cons_vec), name_stats)) in species_kmer_maps
                .iter_mut()
                .zip(species_kmer_consensus.iter_mut())
                .zip(species_kmer_name_stats.iter_mut())
                .enumerate()
                .take(n)
            {
                let present = p.species.get(sid).map(|b| *b).unwrap_or(false);
                if !present {
                    continue;
                }

                for &km in &kmers {
                    match map.entry(km) {
                        Entry::Vacant(e) => {
                            let idx = cons_vec.len();
                            cons_vec.push(None); // consensus will be filled in phase 2
                            name_stats.push(KmerNameStats::default());
                            e.insert(idx);
                        }
                        Entry::Occupied(_) => {
                            // already tracked; nothing to do
                        }
                    }
                }
            }
        }
    }

    // --------------------------------------------------------------------
    // PHASE 2: stream proteomes from disk and build per-k-mer 20AA consensus
    // --------------------------------------------------------------------
    // Get FASTA files and map them to species IDs (same logic as load_proteomes).
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

    let mask = kmer_mask(sym_bits, scan_k);

    struct SpeciesWork {
        sid: usize,
        path: PathBuf,
        map: HashMap<u64, usize>,
        consensus: Vec<Option<Vec<u8>>>,
        name_stats: Vec<KmerNameStats>,
    }

    struct ThreadBuffers {
        aa_ring: Vec<u8>,
        window: Vec<u8>,
        kmer_hits: hashbrown::HashSet<usize>,
    }

    let mut work: Vec<Option<SpeciesWork>> = (0..n).map(|_| None).collect();

    for input in inputs {
        let name = &input.name;
        let Some(&sid) = name_to_sid.get(name) else {
            bail!("Species {:?} from input list not found in graph.", name);
        };

        let map = std::mem::take(&mut species_kmer_maps[sid]);
        let consensus = std::mem::take(&mut species_kmer_consensus[sid]);
        let name_stats = std::mem::take(&mut species_kmer_name_stats[sid]);

        work[sid] = Some(SpeciesWork {
            sid,
            path: input.path.clone(),
            map,
            consensus,
            name_stats,
        });
    }

    let word_interner = Arc::new(Mutex::new(WordInterner::new()));

    let n_threads = num_threads.max(1);
    let pool = ThreadPoolBuilder::new()
        .num_threads(n_threads)
        .build()
        .expect("Failed to build Rayon thread pool in build_species_kmer_consensus");

    pool.install(|| -> Result<()> {
        work.par_iter_mut()
            .filter_map(|opt| opt.as_mut())
            .try_for_each_init(
                || ThreadBuffers {
                    aa_ring: vec![0u8; scan_k],
                    window: Vec::with_capacity(scan_k),
                    kmer_hits: hashbrown::HashSet::new(),
                },
                |buffers, work_item| {
                    if work_item.map.is_empty() {
                        return Ok(());
                    }

                    let reader = open_fasta(&work_item.path)?;
                    let mut fasta = FastaReader::new(reader);

                    while let Some(rec) = fasta.next() {
                        let rec = rec?;
                        let mut name_words: Vec<u32> = Vec::new();
                        let header = match rec.desc() {
                            Some(desc) => {
                                let mut header = rec.id().unwrap_or("").to_string();
                                header.push(' ');
                                header.push_str(desc.unwrap_or(""));
                                header
                            }
                            None => rec.id().unwrap_or("").to_string(),
                        };
                        if !header.is_empty() {
                            let mut interner = word_interner.lock().expect("word interner lock");
                            for (idx, word) in header.split_whitespace().enumerate() {
                                if idx == 0 {
                                    continue;
                                }
                                if word.starts_with('[') {
                                    break;
                                }
                                let word_id = interner.intern(word);
                                if !name_words.contains(&word_id) {
                                    name_words.push(word_id);
                                }
                            }
                        }

                        let seq = rec.seq();
                        if seq.len() < scan_k {
                            continue;
                        }

                        let map = &work_item.map;
                        let cons_vec = &mut work_item.consensus;

                        let mut val: u64 = 0;
                        let mut have: usize = 0;
                        let mut filled: usize = 0;
                        let mut head: usize = 0;

                        for &b in seq {
                            let aa = b.to_ascii_uppercase();
                            let code = recode_byte(aa, recode_scheme);

                            if code == 255 {
                                // Unknown symbol: reset window.
                                val = 0;
                                have = 0;
                                filled = 0;
                                head = 0;
                                continue;
                            }

                            if filled < scan_k {
                                buffers.aa_ring[filled] = aa;
                                filled += 1;
                            } else {
                                buffers.aa_ring[head] = aa;
                                head = (head + 1) % scan_k;
                            }

                            val = (val << sym_bits) | code as u64;

                            if have + 1 < scan_k {
                                have += 1;
                                continue;
                            }

                            if have + 1 > scan_k {
                                val &= mask;
                            } else {
                                have = scan_k;
                            }

                            if have < scan_k {
                                continue;
                            }

                            if let Some(&idx) = map.get(&val) {
                                buffers.window.clear();
                                for i in 0..scan_k {
                                    let pos = (head + i) % scan_k;
                                    buffers.window.push(buffers.aa_ring[pos]);
                                }

                                let entry = &mut cons_vec[idx];
                                if entry.is_none() {
                                    let mut seq_cons = vec![0u8; scan_k];
                                    seq_cons.clone_from_slice(&buffers.window);
                                    *entry = Some(seq_cons);
                                } else {
                                    // Update consensus: first AA wins, conflicts become 'X'.
                                    let seq_cons = entry.as_mut().unwrap();
                                    for (i, &aa_new) in buffers.window.iter().enumerate() {
                                        let aa_old = seq_cons[i];
                                        if aa_old == 0 {
                                            seq_cons[i] = aa_new;
                                        } else if aa_old != aa_new {
                                            seq_cons[i] = b'X';
                                        }
                                    }
                                }
                                if !name_words.is_empty() {
                                    buffers.kmer_hits.insert(idx);
                                }
                            }
                        }

                        if !name_words.is_empty() {
                            for idx in buffers.kmer_hits.drain() {
                                work_item.name_stats[idx].add_words(&name_words);
                            }
                        } else {
                            buffers.kmer_hits.clear();
                        }
                    }

                    // After scanning this species' proteome, turn any still-0 bytes into '-'.
                    for opt in &mut work_item.consensus {
                        if let Some(seq_cons) = opt.as_mut() {
                            for aa in seq_cons.iter_mut() {
                                if *aa == 0 {
                                    *aa = b'-';
                                }
                            }
                        }
                    }

                    Ok(())
                },
            )
    })?;

    for opt in work.into_iter().flatten() {
        species_kmer_maps[opt.sid] = opt.map;
        species_kmer_consensus[opt.sid] = opt.consensus;
        species_kmer_name_stats[opt.sid] = opt.name_stats;
    }

    let word_list = word_interner
        .lock()
        .expect("word interner lock")
        .words
        .clone();

    Ok((
        scan_k,
        species_kmer_maps,
        species_kmer_consensus,
        species_kmer_name_stats,
        word_list,
    ))
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn build_raw_group_block(
    key: (u64, u64),
    paths: &[PathRec],
    k: usize,
    k1: usize,
    head_keep: usize,
    sym_bits: u32,
    scan_k: usize,
    species_kmer_maps: &[HashMap<u64, usize>],
    species_kmer_consensus: &[Vec<Option<Vec<u8>>>],
    species_kmer_name_stats: &[Vec<KmerNameStats>],
    n: usize,
) -> Option<(Vec<Vec<u8>>, KmerNameStats)> {
    let _ = k1;
    let _ = head_keep;
    // Deterministic order for paths: by length, then lexicographically on node IDs
    let mut local_paths: Vec<&PathRec> = paths.iter().collect();
    local_paths.sort_unstable_by(|a, b| {
        a.nodes
            .len()
            .cmp(&b.nodes.len())
            .then_with(|| a.nodes.cmp(&b.nodes))
    });

    // Collect per-path data
    let mut per_path: Vec<PathData> = Vec::with_capacity(local_paths.len());
    let mut block_len: Option<usize> = None;

    for p in local_paths {
        if p.nodes.is_empty() {
            continue;
        }
        let path_codes = render_path_codes(&p.nodes, k, sym_bits);
        if path_codes.is_empty() {
            continue;
        }

        let len = path_codes.len();
        block_len.get_or_insert(len);

        per_path.push(PathData::Recoded {
            codes: path_codes,
            species: &p.species,
        });
    }

    let Some(bl) = block_len else {
        // no usable paths for this group
        return None;
    };

    // One alignment block per group: start with '-' everywhere
    let mut block: Vec<Vec<u8>> = vec![vec![b'-'; bl]; n];

    let mut group_name_stats = KmerNameStats::default();
    let mut group_seen_kmers: hashbrown::HashSet<(usize, usize)> = hashbrown::HashSet::new();

    // Merge paths into the block.
    for data in per_path.into_iter() {
        let PathData::Recoded { codes, species } = data;
        // Repaint from per-k-mer 20AA consensus.
        if scan_k == 0 || codes.len() < scan_k {
            continue;
        }
        let kmers = encoded_kmers_from_codes(&codes, scan_k, sym_bits);

        for (sid, row) in block.iter_mut().enumerate().take(n) {
            if !species.get(sid).map(|b| *b).unwrap_or(false) {
                continue;
            }
            let map = &species_kmer_maps[sid];
            let cons_vec = &species_kmer_consensus[sid];
            let name_vec = &species_kmer_name_stats[sid];

            // Each k-mer kmers[j] covers positions [j .. j+scan_k) in the path.
            for (j, &km) in kmers.iter().enumerate() {
                let Some(&idx) = map.get(&km) else {
                    continue;
                };
                let Some(ref kseq) = cons_vec[idx] else {
                    continue;
                };
                debug_assert_eq!(
                    kseq.len(),
                    scan_k,
                    "Consensus k-mer length mismatch for species {} in group ({}, {})",
                    sid,
                    key.0,
                    key.1
                );
                for (offset, &newc) in kseq.iter().enumerate().take(scan_k) {
                    let pos = j + offset;
                    if pos >= bl {
                        break;
                    }
                    let cur = row[pos];
                    if cur == b'-' {
                        row[pos] = newc;
                    } else if newc == b'-' {
                        // keep cur
                    } else if cur == newc {
                        // identical call; keep cur
                    } else {
                        row[pos] = b'X';
                    }
                }

                if group_seen_kmers.insert((sid, idx)) {
                    let name_stats = &name_vec[idx];
                    if name_stats.total_sets > 0 {
                        group_name_stats.total_sets = group_name_stats
                            .total_sets
                            .saturating_add(name_stats.total_sets);
                        for (word_id, count) in &name_stats.word_counts {
                            if let Some(entry) = group_name_stats
                                .word_counts
                                .iter_mut()
                                .find(|(existing_id, _)| existing_id == word_id)
                            {
                                entry.1 = entry.1.saturating_add(*count);
                            } else {
                                group_name_stats.word_counts.push((*word_id, *count));
                            }
                        }
                    }
                }
            }
        }
    }

    Some((block, group_name_stats))
}
