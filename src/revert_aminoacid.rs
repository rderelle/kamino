use anyhow::Result;
use hashbrown::{HashMap, HashSet};
use seq_io::fasta::{Reader as FastaReader, Record};

use crate::io::{open_fasta, SpeciesInput};
use crate::recode::{recode_byte, RecodeScheme, RECODE_BITS_PER_SYMBOL};
use crate::traverse::{PathRec, VariantGroups};

pub(crate) struct GroupLayout {
    pub key: (u64, u64),
    pub bl: u32,
    pub paths: Vec<GroupPathLayout>,
    pub path_indices: Vec<usize>,
    pub unique_kmers: Vec<u64>,
    pub middle_start: u32,
    pub middle_len: u32,
    pub tail_const_start: u32,
    pub tail_const_keep: u32,
    pub kept_pos: Vec<u16>,
}

pub(crate) struct GroupPathLayout {
    pub kmers: Vec<u64>,
    pub pos0: Vec<u32>,
}

pub(crate) struct SpeciesKmerMap {
    pub map: HashMap<u64, u32>,
    pub windows: Vec<u8>,
    pub name_stats: Vec<KmerNameStats>,
}

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

pub(crate) struct WordInterner {
    map: HashMap<String, u32>,
    pub(crate) words: Vec<String>,
}

impl WordInterner {
    pub(crate) fn new() -> Self {
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

impl GroupLayout {
    #[inline]
    pub(crate) fn kept_index(&self, pos: u32) -> Option<usize> {
        if pos >= self.middle_start && pos < self.middle_start + self.middle_len {
            return Some((pos - self.middle_start) as usize);
        }
        if pos >= self.tail_const_start && pos < self.tail_const_start + self.tail_const_keep {
            let offset = pos - self.tail_const_start;
            return Some((self.middle_len + offset) as usize);
        }
        None
    }

    #[inline]
    pub(crate) fn kept_len(&self) -> usize {
        self.kept_pos.len()
    }
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

pub(crate) fn scan_k_from_k(k: usize) -> usize {
    let min_vg_len: usize = 2 * (k - 1) + 1;
    min_vg_len.min(21)
}
pub(crate) fn build_species_kmer_map(
    input: &SpeciesInput,
    recode_scheme: RecodeScheme,
    scan_k: usize,
    kmers_of_interest: &HashSet<u64>,
    word_interner: &mut WordInterner,
) -> Result<SpeciesKmerMap> {
    if scan_k == 0 {
        return Ok(SpeciesKmerMap {
            map: HashMap::new(),
            windows: Vec::new(),
            name_stats: Vec::new(),
        });
    }

    let sym_bits = RECODE_BITS_PER_SYMBOL;
    let mask = kmer_mask(sym_bits, scan_k);
    let reader = open_fasta(&input.path)?;
    let mut fasta = FastaReader::new(reader);

    let mut map: HashMap<u64, u32> = HashMap::new();
    let mut windows: Vec<u8> = Vec::new();
    let mut name_stats: Vec<KmerNameStats> = Vec::new();

    let mut aa_ring: Vec<u8> = vec![0u8; scan_k];
    let mut window: Vec<u8> = Vec::with_capacity(scan_k);
    let mut kmer_hits: HashSet<u32> = HashSet::new();

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
            for (idx, word) in header.split_whitespace().enumerate() {
                if idx == 0 {
                    continue;
                }
                if word.starts_with('[') {
                    break;
                }
                let word_id = word_interner.intern(word);
                if !name_words.contains(&word_id) {
                    name_words.push(word_id);
                }
            }
        }

        let seq = rec.seq();
        if seq.len() < scan_k {
            continue;
        }

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
                aa_ring[filled] = aa;
                filled += 1;
            } else {
                aa_ring[head] = aa;
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

            if !kmers_of_interest.contains(&val) {
                continue;
            }

            let idx = if let Some(&idx) = map.get(&val) {
                idx
            } else {
                let idx = map.len() as u32;
                map.insert(val, idx);
                windows.resize(windows.len() + scan_k, 0u8);
                name_stats.push(KmerNameStats::default());
                idx
            };

            window.clear();
            for i in 0..scan_k {
                let pos = (head + i) % scan_k;
                window.push(aa_ring[pos]);
            }

            let start = idx as usize * scan_k;
            let slot = &mut windows[start..start + scan_k];
            for (i, &aa_new) in window.iter().enumerate() {
                let aa_old = slot[i];
                if aa_old == 0 {
                    slot[i] = aa_new;
                } else if aa_old != aa_new {
                    slot[i] = b'X';
                }
            }
            if !name_words.is_empty() {
                kmer_hits.insert(idx);
            }
        }

        if !name_words.is_empty() {
            for idx in kmer_hits.drain() {
                name_stats[idx as usize].add_words(&name_words);
            }
        } else {
            kmer_hits.clear();
        }
    }

    for slot in windows.iter_mut() {
        if *slot == 0 {
            *slot = b'-';
        }
    }

    Ok(SpeciesKmerMap {
        map,
        windows,
        name_stats,
    })
}

pub(crate) fn build_group_layouts(
    groups: &VariantGroups,
    k: usize,
    head_keep: usize,
    scan_k: usize,
) -> Vec<GroupLayout> {
    let sym_bits = RECODE_BITS_PER_SYMBOL;
    let k1 = k - 1;

    let mut keys: Vec<(u64, u64)> = groups.keys().copied().collect();
    keys.sort_unstable();

    let mut layouts = Vec::with_capacity(keys.len());

    for key in keys {
        let paths = &groups[&key];
        let mut local_paths: Vec<(usize, &PathRec)> = paths.iter().enumerate().collect();
        local_paths.sort_unstable_by(|(_, a), (_, b)| {
            a.nodes
                .len()
                .cmp(&b.nodes.len())
                .then_with(|| a.nodes.cmp(&b.nodes))
        });

        let mut paths_layout: Vec<GroupPathLayout> = Vec::with_capacity(local_paths.len());
        let mut path_indices: Vec<usize> = Vec::with_capacity(local_paths.len());
        let mut block_len: Option<usize> = None;
        let mut unique_kmers: HashSet<u64> = HashSet::new();

        for (idx, p) in local_paths {
            let mut path_kmers_out: Vec<u64> = Vec::new();
            let mut path_pos0_out: Vec<u32> = Vec::new();
            if p.nodes.is_empty() {
                continue;
            }
            let path_codes = render_path_codes(&p.nodes, k, sym_bits);
            if path_codes.is_empty() {
                continue;
            }
            let len = path_codes.len();
            block_len.get_or_insert(len);

            if scan_k == 0 || path_codes.len() < scan_k {
                paths_layout.push(GroupPathLayout {
                    kmers: path_kmers_out,
                    pos0: path_pos0_out,
                });
                path_indices.push(idx);
                continue;
            }
            let path_kmers = encoded_kmers_from_codes(&path_codes, scan_k, sym_bits);
            for (j, &km) in path_kmers.iter().enumerate() {
                path_kmers_out.push(km);
                path_pos0_out.push(j as u32);
                unique_kmers.insert(km);
            }

            paths_layout.push(GroupPathLayout {
                kmers: path_kmers_out,
                pos0: path_pos0_out,
            });
            path_indices.push(idx);
        }

        let Some(bl) = block_len else {
            continue;
        };

        let middle_start = k1 as u32;
        let middle_end = bl.saturating_sub(k1) as u32;
        let middle_len = middle_end.saturating_sub(middle_start);
        let tail_const_start = middle_end;
        let tail_const_keep = head_keep.min(bl.saturating_sub(k1)) as u32;

        let mut kept_pos: Vec<u16> =
            Vec::with_capacity((middle_len + tail_const_keep).try_into().unwrap_or(0usize));
        for pos in middle_start..middle_end {
            kept_pos.push(pos as u16);
        }
        for pos in tail_const_start..tail_const_start + tail_const_keep {
            kept_pos.push(pos as u16);
        }

        layouts.push(GroupLayout {
            key,
            bl: bl as u32,
            paths: paths_layout,
            path_indices,
            unique_kmers: unique_kmers.into_iter().collect(),
            middle_start,
            middle_len,
            tail_const_start,
            tail_const_keep,
            kept_pos,
        });
    }

    layouts
}
