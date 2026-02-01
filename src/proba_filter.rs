//! Probabilistic prefilter for proteins using:
//! - per-proteome Bloom filter to count each (k-1)-mer once per species.
//! - global Count-Min Sketch (CMS) to estimate how many species contain each (k-1)-mer.
//!
//! False positives tradeoffs:
//! - Bloom false positives cause undercounting ((k-1)-mers may look already seen),
//!   which can drop some proteins that should pass (conservative).
//! - CMS overestimates counts, which can let extra proteins through (permissive).
//!   The combination keeps sensitivity while bounding memory and runtime.

use crate::recode::{recode_byte, RecodeScheme};
use std::sync::atomic::{AtomicU16, Ordering};

pub const BLOOM_BITS: usize = 1 << 23;
pub const BLOOM_HASHES: usize = 2;
pub const CMS_WIDTH: usize = 1 << 24;
pub const CMS_DEPTH: usize = 4;
pub const MIN_GOOD_KMERS: u32 = 2;

const SEEDS: [u64; 6] = [
    0x9e3779b97f4a7c15,
    0xbf58476d1ce4e5b9,
    0x94d049bb133111eb,
    0x3c79ac492ba7b653,
    0x1c69b3f74ac4ae35,
    0x6a09e667f3bcc909,
];

fn mix64(mut x: u64) -> u64 {
    x ^= x >> 30;
    x = x.wrapping_mul(0xbf58476d1ce4e5b9);
    x ^= x >> 27;
    x = x.wrapping_mul(0x94d049bb133111eb);
    x ^ (x >> 31)
}

#[derive(Clone)]
pub struct BloomFilter {
    bits: Vec<u64>,
    mask: u64,
}

impl BloomFilter {
    pub fn new(bits: usize) -> Self {
        let size = bits.max(64).next_power_of_two();
        let words = size.div_ceil(64);
        BloomFilter {
            bits: vec![0; words],
            mask: (size - 1) as u64,
        }
    }

    /// Returns true when the item is likely already present.
    pub fn test_and_set(&mut self, key: u64) -> bool {
        let h1 = mix64(key ^ SEEDS[0]);
        let h2 = mix64(key ^ SEEDS[1]) | 1;
        let mut seen = true;
        for i in 0..BLOOM_HASHES {
            let idx = (h1.wrapping_add((i as u64).wrapping_mul(h2)) & self.mask) as usize;
            let word = idx >> 6;
            let bit = 1u64 << (idx & 63);
            let had = (self.bits[word] & bit) != 0;
            self.bits[word] |= bit;
            if !had {
                seen = false;
            }
        }
        seen
    }
}

#[derive(Clone)]
pub struct CountMinSketch {
    width: usize,
    table: Vec<u16>,
    seeds: Vec<u64>,
}

impl CountMinSketch {
    pub fn new(width: usize, depth: usize) -> Self {
        let w = width.max(1).next_power_of_two();
        let d = depth.max(1);
        let mut seeds = Vec::with_capacity(d);
        for i in 0..d {
            seeds.push(SEEDS[2 + (i % (SEEDS.len() - 2))]);
        }
        CountMinSketch {
            width: w,
            table: vec![0; w * d],
            seeds,
        }
    }

    pub fn estimate(&self, key: u64) -> u32 {
        let mut min = u16::MAX;
        for (row, seed) in self.seeds.iter().enumerate() {
            let idx = (mix64(key ^ seed) & (self.width as u64 - 1)) as usize;
            let offset = row * self.width + idx;
            min = min.min(self.table[offset]);
        }
        min as u32
    }
}

/// Thread-safe CMS used during parallel ingestion.
///
/// Each (k-1)-mer is incremented at most once per species by guarding with a
/// per-species Bloom filter, so the counts represent species occurrences.
pub struct AtomicCountMinSketch {
    width: usize,
    table: Vec<AtomicU16>,
    seeds: Vec<u64>,
}

impl AtomicCountMinSketch {
    pub fn new(width: usize, depth: usize) -> Self {
        let w = width.max(1).next_power_of_two();
        let d = depth.max(1);
        let mut seeds = Vec::with_capacity(d);
        for i in 0..d {
            seeds.push(SEEDS[2 + (i % (SEEDS.len() - 2))]);
        }
        let mut table = Vec::with_capacity(w * d);
        table.resize_with(w * d, || AtomicU16::new(0));
        AtomicCountMinSketch {
            width: w,
            table,
            seeds,
        }
    }

    /// Increment all CMS rows for a key using saturating arithmetic.
    pub fn increment(&self, key: u64) {
        for (row, seed) in self.seeds.iter().enumerate() {
            let idx = (mix64(key ^ seed) & (self.width as u64 - 1)) as usize;
            let offset = row * self.width + idx;
            let cell = &self.table[offset];
            let mut current = cell.load(Ordering::Relaxed);
            loop {
                if current == u16::MAX {
                    break;
                }
                let next = current.saturating_add(1);
                match cell.compare_exchange_weak(
                    current,
                    next,
                    Ordering::Relaxed,
                    Ordering::Relaxed,
                ) {
                    Ok(_) => break,
                    Err(actual) => current = actual,
                }
            }
        }
    }

    /// Create a read-only snapshot for the filtering pass.
    pub fn snapshot(&self) -> CountMinSketch {
        let mut cms = CountMinSketch::new(self.width, self.seeds.len());
        for (slot, cell) in cms.table.iter_mut().zip(self.table.iter()) {
            *slot = cell.load(Ordering::Relaxed);
        }
        cms
    }
}

pub fn stream_recoded_kmers(
    seq: &[u8],
    k: usize,
    sym_bits: u32,
    k_mask: u64,
    scheme: RecodeScheme,
    mut on_kmer: impl FnMut(u64),
) {
    let mut roll: u64 = 0;
    let mut have: usize = 0;
    for &b in seq {
        let a = recode_byte(b, scheme);
        if a == 255 {
            roll = 0;
            have = 0;
            continue;
        }
        roll = ((roll << sym_bits) | (a as u64)) & k_mask;
        if have < k {
            have += 1;
        }
        if have == k {
            on_kmer(roll);
        }
    }
}

#[allow(clippy::too_many_arguments)]
pub fn stream_recoded_edges(
    seq: &[u8],
    k1: usize,
    sym_bits: u32,
    k1_mask: u64,
    scheme: RecodeScheme,
    mut on_edge: impl FnMut(u64, u32, u64),
) {
    let mut roll: u64 = 0;
    let mut have: usize = 0;
    let mut prev: Option<u64> = None;
    for &b in seq {
        let a = recode_byte(b, scheme);
        if a == 255 {
            roll = 0;
            have = 0;
            prev = None;
            continue;
        }
        roll = ((roll << sym_bits) | (a as u64)) & k1_mask;
        if have < k1 {
            have += 1;
        }
        if have == k1 {
            if let Some(u) = prev {
                let sym = a as u32;
                on_edge(u, sym, roll);
            }
            prev = Some(roll);
        }
    }
}

#[allow(clippy::too_many_arguments)]
pub fn protein_passes_counts(
    seq: &[u8],
    k: usize,
    sym_bits: u32,
    k_mask: u64,
    scheme: RecodeScheme,
    min_needed: u32,
    length_middle: usize,
    mut count: impl FnMut(u64) -> u32,
) -> bool {
    let mut roll: u64 = 0;
    let mut have: usize = 0;
    let mut last_hit_start: Option<usize> = None;
    let mut good_hits = 0u32;
    for (idx, &b) in seq.iter().enumerate() {
        let a = recode_byte(b, scheme);
        if a == 255 {
            roll = 0;
            have = 0;
            continue;
        }
        roll = ((roll << sym_bits) | (a as u64)) & k_mask;
        if have < k {
            have += 1;
        }
        if have == k && count(roll) >= min_needed {
            let start = idx + 1 - k;
            if let Some(prev_start) = last_hit_start {
                if start >= prev_start + k {
                    let gap = start - (prev_start + k);
                    if gap <= length_middle {
                        good_hits += 1;
                        if good_hits >= MIN_GOOD_KMERS {
                            return true;
                        }
                    } else {
                        good_hits = 1;
                    }
                    last_hit_start = Some(start);
                }
            } else {
                good_hits = 1;
                last_hit_start = Some(start);
            }
        }
    }
    false
}
