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

#[inline]
fn four_offsets(key: u64, width: usize) -> [usize; 4] {
    let mask = width as u64 - 1;
    [
        (mix64(key ^ SEEDS[2]) & mask) as usize,
        width + (mix64(key ^ SEEDS[3]) & mask) as usize,
        2 * width + (mix64(key ^ SEEDS[4]) & mask) as usize,
        3 * width + (mix64(key ^ SEEDS[5]) & mask) as usize,
    ]
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

    /// Single-writer counting avoids atomic read-modify-write operations.
    pub fn increment(&mut self, key: u64) {
        for (row, seed) in self.seeds.iter().enumerate() {
            let idx = (mix64(key ^ seed) & (self.width as u64 - 1)) as usize;
            let cell = &mut self.table[row * self.width + idx];
            *cell = cell.saturating_add(1);
        }
    }

    /// Return the exact occupancy for shared anchors. For an unshared k-mer,
    /// any value below `min_needed` preserves all neighbour-boundary comparisons,
    /// so the lookup stops at the first row below the threshold.
    pub fn estimate_for_anchors(&self, key: u64, min_needed: u32) -> u32 {
        // The production sketch has four rows: computing the four probes up front
        // lets their loads proceed independently instead of looping over `seeds`.
        if self.seeds.len() == 4 {
            let [a, b, c, d] = four_offsets(key, self.width);
            let first = self.table[a];
            if (first as u32) < min_needed {
                return first as u32;
            }
            return first
                .min(self.table[b])
                .min(self.table[c])
                .min(self.table[d]) as u32;
        }
        let mut min = u16::MAX;
        for (row, seed) in self.seeds.iter().enumerate() {
            let idx = (mix64(key ^ seed) & (self.width as u64 - 1)) as usize;
            let value = self.table[row * self.width + idx];
            if (value as u32) < min_needed {
                return value as u32;
            }
            min = min.min(value);
        }
        min as u32
    }

    #[cfg(test)]
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

    /// Consume the finished sketch without retaining a second 512 MiB table.
    /// Owned iteration allows Vec to reuse the allocation when converting counters.
    pub fn into_snapshot(self) -> CountMinSketch {
        CountMinSketch {
            width: self.width,
            table: self.table.into_iter().map(AtomicU16::into_inner).collect(),
            seeds: self.seeds,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn plain_and_atomic_counts_match_including_collisions_and_saturation() {
        for width in [1, 16, 1024] {
            for depth in [1, 2, 4, 7] {
                let mut plain = CountMinSketch::new(width, depth);
                let atomic = AtomicCountMinSketch::new(width, depth);
                for key in (0..1000).map(|i| mix64(i % 311)) {
                    plain.increment(key);
                    atomic.increment(key);
                }
                // Exercise the saturation boundary without a large input file.
                for _ in 0..65536 {
                    plain.increment(42);
                    atomic.increment(42);
                }
                let frozen = atomic.into_snapshot();
                assert_eq!(plain.table, frozen.table);
                for key in (0..1200).map(mix64) {
                    let reference = plain
                        .seeds
                        .iter()
                        .enumerate()
                        .map(|(row, seed)| {
                            plain.table[row * plain.width
                                + (mix64(key ^ seed) & (plain.width as u64 - 1)) as usize]
                        })
                        .min()
                        .unwrap() as u32;
                    assert_eq!(plain.estimate(key), reference);
                    assert_eq!(frozen.estimate(key), reference);
                    for threshold in [1, 2, 3, 10, 65535] {
                        let observed = frozen.estimate_for_anchors(key, threshold);
                        assert_eq!(observed >= threshold, reference >= threshold);
                        if reference >= threshold {
                            assert_eq!(observed, reference);
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn below_threshold_counts_do_not_change_anchor_boundaries() {
        for threshold in 1..=32 {
            for left in 0..=32 {
                for right in 0..=32 {
                    let boundaries = |a, b| (a >= threshold && a > b, b >= threshold && a < b);
                    let clamp = |value| {
                        if value < threshold {
                            threshold - 1
                        } else {
                            value
                        }
                    };
                    assert_eq!(
                        boundaries(left, right),
                        boundaries(clamp(left), clamp(right))
                    );
                }
            }
        }
    }

    #[test]
    fn concurrent_updates_keep_saturating_counts() {
        let atomic = AtomicCountMinSketch::new(1, 4);
        std::thread::scope(|scope| {
            for _ in 0..4 {
                let atomic = &atomic;
                scope.spawn(move || {
                    for _ in 0..20000 {
                        atomic.increment(42);
                    }
                });
            }
        });
        assert_eq!(atomic.into_snapshot().estimate(42), u16::MAX as u32);
    }
}
