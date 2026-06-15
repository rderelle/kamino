//! Probabilistic prefilter for proteins using:
//! - per-proteome Bloom filter to count each k-mer once per species.
//! - global Count-Min Sketch (CMS) to estimate how many species contain each k-mer.
//!
//! False positives tradeoffs:
//! - Bloom false positives cause undercounting (k-mers may look already seen),
//!   which can drop some proteins that should pass (conservative).
//! - CMS overestimates counts, which can let extra proteins through (permissive).
//!   The combination keeps sensitivity while bounding memory and runtime.

use std::sync::atomic::{AtomicU16, Ordering};

pub const BLOOM_BITS: usize = 1 << 25;
/// Number of independent probes used by each per-species Bloom filter.
pub const BLOOM_HASHES: usize = 2;
/// Number of counters per count-min-sketch row; must stay a power-of-two size.
pub const CMS_WIDTH: usize = 1 << 26;
/// Number of hash rows in the count-min sketch.
pub const CMS_DEPTH: usize = 4;

const SEEDS: [u64; 6] = [
    0x9e3779b97f4a7c15,
    0xbf58476d1ce4e5b9,
    0x94d049bb133111eb,
    0x3c79ac492ba7b653,
    0x1c69b3f74ac4ae35,
    0x6a09e667f3bcc909,
];

fn mix64(mut x: u64) -> u64 {
    // SplitMix64 finalizer: cheap avalanche hash for already-encoded k-mer keys.
    x ^= x >> 30;
    x = x.wrapping_mul(0xbf58476d1ce4e5b9);
    x ^= x >> 27;
    x = x.wrapping_mul(0x94d049bb133111eb);
    x ^ (x >> 31)
}

#[derive(Clone)]
pub struct BloomFilter {
    /// Packed bitset for the filter.
    bits: Vec<u64>,
    /// Power-of-two mask used instead of modulo when mapping hashes to bit indexes.
    mask: u64,
}

impl BloomFilter {
    pub fn new(bits: usize) -> Self {
        // Round up to a power of two so indexing can use `hash & mask`.
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
    /// Counters per row.
    width: usize,
    /// Row-major counter table.
    table: Vec<u16>,
    /// Per-row seeds that produce independent hash locations.
    seeds: Vec<u64>,
}

impl CountMinSketch {
    pub fn new(width: usize, depth: usize) -> Self {
        // Use at least one row/counter even if a caller accidentally passes zero.
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
        // CMS query returns the minimum row counter, which bounds overestimation.
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
/// Each k-mer is incremented at most once per species by guarding with a
/// per-species Bloom filter, so the counts represent species occurrences.
pub struct AtomicCountMinSketch {
    /// Counters per row.
    width: usize,
    /// Atomic row-major table for parallel counting.
    table: Vec<AtomicU16>,
    /// Per-row seeds matching the read-only snapshot representation.
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
