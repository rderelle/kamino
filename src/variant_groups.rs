use anyhow::{bail, Result};
use bitvec::prelude::*;
use hashbrown::HashMap;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use seq_io::fasta::{Reader as FastaReader, Record};
use std::path::{Path, PathBuf};

use crate::graph::Graph;
use crate::io::{collect_fasta_files, open_fasta, species_name_from_path};
use crate::recode::recode_byte;
use crate::traverse::{PathRec, VariantGroups};

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
    pub dropped_blocks_due_to_middle_filter: usize,
    pub dropped_blocks_due_to_length: usize,
}

/// Outcome of attempting to build a single variant-group alignment block.
pub(crate) enum GroupOutcome {
    Kept {
        filtered_rows: Vec<Vec<u8>>,
        filtered_len: usize,
    },
    DroppedMiddle,
    DroppedLength,
    Skipped,
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
pub(crate) fn build_group_block(
    key: (u64, u64),
    paths: &[PathRec],
    k: usize,
    k1: usize,
    head_keep: usize,
    sym_bits: u32,
    scan_k: usize,
    species_kmer_maps: &[HashMap<u64, usize>],
    species_kmer_consensus: &[Vec<Option<Vec<u8>>>],
    bubble_ratio: f64,
    max_middle_len: usize,
    n: usize,
) -> GroupOutcome {
    // Produce a single alignment block for the group identified by `key`,
    // repainting sequences from per-species consensus and respecting the
    // middle-length and middle-uniqueness filters.
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
        return GroupOutcome::Skipped;
    };

    // One alignment block per group: start with '-' everywhere
    let mut block: Vec<Vec<u8>> = vec![vec![b'-'; bl]; n];

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
            }
        }
    }

    // Region indices
    let lead_full = k1;
    let trail_full = k1;

    let middle_start = lead_full;
    let middle_end = bl - trail_full; // exclusive
    let middle_len = middle_end - middle_start;

    // Missing cutoff based on presence ratio requirement (presence_ratio = 1 - missing_ratio >= bubble_ratio)
    let missing_cutoff = 1.0 - bubble_ratio;

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
        return GroupOutcome::DroppedMiddle;
    }
    
    // If no middle column is polymorphic, drop the block.
    if !has_polymorphic_middle {
        return GroupOutcome::DroppedMiddle;
    }
    
    if keep_middle_cols.len() > max_middle_len {
        return GroupOutcome::DroppedLength;
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

    GroupOutcome::Kept {
        filtered_rows,
        filtered_len: block_len_filtered,
    }
}

#[allow(clippy::too_many_arguments)]
pub fn build_variant_group_alignment(
    input_dir: &Path,
    g: &Graph,
    groups: &VariantGroups,
    k: usize,
    head_max: usize,
    bubble_ratio: f32,
    max_middle_len: usize,
    num_threads: usize,
) -> Result<VariantGroupAlignment> {
    let species = &g.species_names;
    let n = g.n_species;

    let k1 = k - 1;
    let head_keep = head_max;
    let sym_bits = g.sym_bits;

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

    for paths in groups.values() {
        for p in paths {
            // Recoded symbol sequence for this path (Dayhoff6 codes).
            let path_codes = render_path_codes(&p.nodes, k, sym_bits);
            let kmers = encoded_kmers_from_codes(&path_codes, scan_k, sym_bits);

            // For every species present on this path, register these k-mers.
            for (sid, (map, cons_vec)) in species_kmer_maps
                .iter_mut()
                .zip(species_kmer_consensus.iter_mut())
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
    let files = collect_fasta_files(input_dir)?;
    if files.len() != species.len() {
        bail!(
            "Input FASTA count ({}) does not match graph species count ({}).",
            files.len(),
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
    }

    struct ThreadBuffers {
        aa_ring: Vec<u8>,
        window: Vec<u8>,
    }

    let mut work: Vec<Option<SpeciesWork>> = (0..n).map(|_| None).collect();

    for path in files {
        let name = species_name_from_path(&path);
        let Some(&sid) = name_to_sid.get(&name) else {
            bail!(
                "Species {:?} from input directory not found in graph.",
                name
            );
        };

        let map = std::mem::take(&mut species_kmer_maps[sid]);
        let consensus = std::mem::take(&mut species_kmer_consensus[sid]);

        work[sid] = Some(SpeciesWork {
            sid,
            path,
            map,
            consensus,
        });
    }

    let n_threads = num_threads.max(1);
    let pool = ThreadPoolBuilder::new()
        .num_threads(n_threads)
        .build()
        .expect("Failed to build Rayon thread pool in build_variant_group_alignment");

    pool.install(|| -> Result<()> {
        work.par_iter_mut()
            .filter_map(|opt| opt.as_mut())
            .try_for_each_init(
                || ThreadBuffers {
                    aa_ring: vec![0u8; scan_k],
                    window: Vec::with_capacity(scan_k),
                },
                |buffers, work_item| {
                    if work_item.map.is_empty() {
                        return Ok(());
                    }

                    let reader = open_fasta(&work_item.path)?;
                    let mut fasta = FastaReader::new(reader);

                    while let Some(rec) = fasta.next() {
                        let rec = rec?;
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
                            let code = recode_byte(aa);

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
                            }
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
    }

    crate::output::build_concatenated_alignment(
        groups,
        k,
        k1,
        head_keep,
        sym_bits,
        scan_k,
        &species_kmer_maps,
        &species_kmer_consensus,
        bubble_ratio,
        max_middle_len,
        n,
    )
}
