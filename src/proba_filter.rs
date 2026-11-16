//! Probabilistic prefilter for proteins using:
//! - per-proteome Bloom filter to count each k-mer once per species.
//! - global Count-Min Sketch (CMS) to estimate how many species contain each k-mer.
//!
//! False positives tradeoffs:
//! - Bloom false positives cause undercounting (k-mers may look already seen),
//!   which can drop some proteins that should pass (conservative).
//! - CMS overestimates counts, which can let extra proteins through (permissive).
//!   The combination keeps sensitivity while bounding memory and runtime.
//!
//! Both tables are far larger than the CPU caches, and each k-mer reads them at
//! unrelated places. K-mers are therefore handled in batches, and the memory needed
//! by a k-mer is prefetched [`AHEAD`] k-mers before it is used: otherwise most of
//! the time is spent waiting for memory, one k-mer at a time.

use std::sync::atomic::{AtomicU16, Ordering};

pub const BLOOM_BITS: usize = 1 << 25;
/// Number of independent probes used by each per-species Bloom filter.
pub const BLOOM_HASHES: usize = 2;
/// Number of counters per count-min-sketch row; must stay a power-of-two size.
pub const CMS_WIDTH: usize = 1 << 26;
/// Number of hash rows in the count-min sketch.
pub const CMS_DEPTH: usize = 4;

/// K-mers handed to the Bloom filter and to the sketch in one call.
pub const BATCH: usize = 1 << 12;
/// Number of k-mers between prefetching the memory of a k-mer and using it.
const AHEAD: usize = 32;

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

/// Ask the CPU to start loading the cache line of `target`. This is only a hint: it
/// never changes results, and does nothing on other architectures.
#[inline(always)]
fn prefetch<T>(target: &T) {
    let p: *const u8 = (target as *const T).cast();
    // SAFETY (both blocks): a prefetch is not a memory access for the program (it
    // cannot fault), and `p` comes from a valid reference anyway.
    #[cfg(target_arch = "x86_64")]
    unsafe {
        use std::arch::x86_64::{_mm_prefetch, _MM_HINT_T0};
        _mm_prefetch(p.cast(), _MM_HINT_T0);
    }
    #[cfg(target_arch = "aarch64")]
    unsafe {
        std::arch::asm!(
            "prfm pldl1keep, [{p}]",
            p = in(reg) p,
            options(nostack, preserves_flags, readonly)
        );
    }
    #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
    let _ = p;
}

/// Table offset of the counter of `key` in each row of a sketch.
#[inline]
fn row_offsets(seeds: &[u64], width: usize, key: u64) -> impl Iterator<Item = usize> + '_ {
    let mask = width as u64 - 1;
    seeds
        .iter()
        .enumerate()
        .map(move |(row, seed)| row * width + (mix64(key ^ seed) & mask) as usize)
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

    /// Bit indexes probed for a key.
    #[inline]
    fn probes(&self, key: u64) -> [usize; BLOOM_HASHES] {
        let h1 = mix64(key ^ SEEDS[0]);
        let h2 = mix64(key ^ SEEDS[1]) | 1;
        std::array::from_fn(|i| (h1.wrapping_add((i as u64).wrapping_mul(h2)) & self.mask) as usize)
    }

    /// Returns true when the item is likely already present.
    pub fn test_and_set(&mut self, key: u64) -> bool {
        let mut seen = true;
        for idx in self.probes(key) {
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

    /// Keep, in order, the keys that were not present yet: the same decisions as
    /// [`Self::test_and_set`] on each key in turn.
    pub fn retain_unseen(&mut self, keys: &mut Vec<u64>) {
        let mut kept = 0;
        for i in 0..keys.len() {
            if let Some(&later) = keys.get(i + AHEAD) {
                for idx in self.probes(later) {
                    prefetch(&self.bits[idx >> 6]);
                }
            }
            let key = keys[i];
            if !self.test_and_set(key) {
                keys[kept] = key;
                kept += 1;
            }
        }
        keys.truncate(kept);
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

    /// [`Self::increment`] for each key, its counters being prefetched beforehand.
    pub fn increment_batch(&mut self, keys: &[u64]) {
        for (i, &key) in keys.iter().enumerate() {
            if let Some(&later) = keys.get(i + AHEAD) {
                row_offsets(&self.seeds, self.width, later).for_each(|o| prefetch(&self.table[o]));
            }
            self.increment(key);
        }
    }

    /// [`Self::estimate_for_anchors`] for each key, in order.
    ///
    /// The first-row counter of a key is prefetched [`AHEAD`] keys in advance. Halfway
    /// there, it has arrived and is read: the other counters of the key are needed,
    /// and thus prefetched in turn, only if this first one reaches `min_needed`.
    pub fn estimate_batch_for_anchors(&self, keys: &[u64], min_needed: u32) -> Vec<u32> {
        let first_row = |key| (mix64(key ^ self.seeds[0]) & (self.width as u64 - 1)) as usize;
        let mut estimates = Vec::with_capacity(keys.len());
        for (i, &key) in keys.iter().enumerate() {
            if let Some(&later) = keys.get(i + AHEAD) {
                prefetch(&self.table[first_row(later)]);
            }
            if let Some(&soon) = keys.get(i + AHEAD / 2) {
                if self.table[first_row(soon)] as u32 >= min_needed {
                    row_offsets(&self.seeds, self.width, soon)
                        .skip(1)
                        .for_each(|o| prefetch(&self.table[o]));
                }
            }
            estimates.push(self.estimate_for_anchors(key, min_needed));
        }
        estimates
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

    /// [`Self::increment`] for each key, its counters being prefetched beforehand.
    pub fn increment_batch(&self, keys: &[u64]) {
        for (i, &key) in keys.iter().enumerate() {
            if let Some(&later) = keys.get(i + AHEAD) {
                row_offsets(&self.seeds, self.width, later).for_each(|o| prefetch(&self.table[o]));
            }
            self.increment(key);
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
    fn batched_bloom_makes_the_same_decisions_as_single_keys() {
        // Small filter and repeated keys: many keys are (or look) already present.
        let keys: Vec<u64> = (0..5000u64).map(|i| mix64(i % 1777)).collect();
        let mut single = BloomFilter::new(1 << 12);
        let expected: Vec<u64> = keys
            .iter()
            .copied()
            .filter(|&key| !single.test_and_set(key))
            .collect();
        assert!(expected.len() < 1777);
        let mut batched = BloomFilter::new(1 << 12);
        let mut kept = Vec::new();
        for chunk in keys.chunks(700) {
            let mut chunk = chunk.to_vec();
            batched.retain_unseen(&mut chunk);
            kept.extend(chunk);
        }
        assert_eq!(kept, expected);
        assert_eq!(batched.bits, single.bits);
    }

    #[test]
    fn batched_sketch_calls_match_single_key_calls() {
        let keys: Vec<u64> = (0..3000u64).map(|i| mix64(i % 701)).collect();
        for width in [1, 16, 1024] {
            for depth in [1, 2, 4, 7] {
                let mut single = CountMinSketch::new(width, depth);
                keys.iter().for_each(|&key| single.increment(key));
                let mut plain = CountMinSketch::new(width, depth);
                let atomic = AtomicCountMinSketch::new(width, depth);
                for chunk in keys.chunks(1000) {
                    plain.increment_batch(chunk);
                    atomic.increment_batch(chunk);
                }
                assert_eq!(plain.table, single.table);
                assert_eq!(atomic.into_snapshot().table, single.table);
                // Counts are 4 or 5 per key, plus collisions.
                let queries: Vec<u64> = (0..1500).map(mix64).collect();
                for threshold in [1, 4, 5, 6, 9, 10, 65535] {
                    let expected: Vec<u32> = queries
                        .iter()
                        .map(|&key| single.estimate_for_anchors(key, threshold))
                        .collect();
                    assert_eq!(
                        single.estimate_batch_for_anchors(&queries, threshold),
                        expected
                    );
                }
                assert!(single.estimate_batch_for_anchors(&[], 1).is_empty());
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
