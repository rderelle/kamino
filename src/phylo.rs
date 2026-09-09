//! Minimal neighbor-joining implementation used for the optional `--nj` tree.
//!
//! Distances are computed from pairwise amino-acid differences while ignoring gaps
//! and unknown residues; the resulting matrix is converted to Newick with the
//! classic neighbor-joining reduction loop. Optional non-parametric bootstrap
//! support is calculated by resampling alignment columns with replacement.
//!
//! Four things keep this fast enough for thousands of taxa and a million columns:
//!
//! * The alignment is **bit-sliced**: six u64 planes per taxon per 64 columns (five
//!   residue-code bits plus a validity bit), so a taxon pair is compared 64 columns
//!   at a time with a handful of XOR/AND/`count_ones` instead of a branch per column.
//! * Columns that are constant *and* gap-free are a valid match for *every* pair, and
//!   all-gap columns affect no pair. Both are folded into one scalar per replicate and
//!   dropped from the per-pair loops. This is exact, not an approximation.
//! * Bootstrap column multiplicities are stored as **bit planes**, so a replicate's
//!   weighted site counts are `Σ 2^p · popcount(mask & plane_p)`: about `4·L/64` word
//!   operations per pair instead of one step per sampled column, and ~10x less memory.
//! * The reduction runs on a compact n×n matrix with swap-removal (both scans read
//!   memory sequentially), maintains row sums incrementally, and is executed **in
//!   parallel across replicates**.
//!
//! Output is deterministic and independent of the thread count: exact ties in the Q
//! criterion are broken on taxon index rather than on iteration order.
use rayon::prelude::*;
use rayon::{ThreadPool, ThreadPoolBuilder};
use std::collections::HashMap;
use std::fmt::Write as _;
use std::sync::Mutex;

#[derive(Debug, Clone)]
pub struct Edge {
    /// Child node index in the shared node arena.
    pub child: usize,
    /// Branch length from parent to child.
    pub len: f64,
}

#[derive(Debug, Clone)]
pub struct Node {
    /// Leaf name; internal nodes use `None`.
    pub name: Option<String>,
    /// Outgoing edges from this rooted representation.
    pub children: Vec<Edge>,
}

// LG stationary frequencies from Le & Gascuel (2008).
const LG_PI: [f64; 20] = [
    0.079_066, 0.055_941, 0.041_977, 0.053_052, 0.012_937, 0.040_767, 0.071_586, 0.057_337,
    0.022_355, 0.062_157, 0.099_081, 0.064_600, 0.022_951, 0.042_302, 0.044_040, 0.061_197,
    0.053_287, 0.012_777, 0.027_843, 0.070_200,
];

// Pre-computed constant c = 1.0 - Σπᵢ² (hoisted out of the hot loop).
const C: f64 = {
    let mut s = 0.0f64;
    let mut i = 0usize;
    while i < LG_PI.len() {
        let pi = LG_PI[i];
        s += pi * pi;
        i += 1;
    }
    1.0 - s
};

// Fixed seed keeps bootstrap output reproducible across runs without adding another CLI option.
const BOOTSTRAP_SEED: u64 = 0x6b61_6d69_6e6f_0001;

/// Small deterministic PRNG (SplitMix64) used only for bootstrap column sampling.
struct BootstrapRng {
    state: u64,
}

impl BootstrapRng {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }

    fn sample_index(&mut self, upper: usize) -> usize {
        debug_assert!(upper > 0);
        let upper = upper as u64;
        // Rejection sampling avoids modulo bias.
        let zone = u64::MAX - (u64::MAX % upper);
        loop {
            let x = self.next_u64();
            if x < zone {
                return (x % upper) as usize;
            }
        }
    }
}

fn corrected_distance(valid_sites: u64, mismatches: u64) -> f64 {
    if valid_sites == 0 {
        return 0.0;
    }

    let p = (mismatches as f64) / (valid_sites as f64);
    let mut pc = p / C;
    if pc >= 1.0 {
        pc = 0.999_999_999;
    }
    -(1.0 - pc).ln()
}

// ------------------------------
// Bit-sliced, column-reduced alignment
// ------------------------------

/// Code for gaps and unknown residues; anything below it is a real residue.
const GAP: u8 = 20;
/// u64 planes stored per taxon per 64-column block: five code bits plus validity.
const PLANES: usize = 6;
/// Index of the validity plane within a block.
const VALID: usize = 5;
/// `col_map` sentinel: all-gap column, affects no pair.
const DROP: u32 = u32::MAX;
/// `col_map` sentinel: constant gap-free column, a valid match for every pair.
const CONST: u32 = u32::MAX - 1;

/// Byte -> LG residue index (A,R,N,...,V); everything else maps to [`GAP`].
const AA_CODE: [u8; 256] = {
    let mut table = [GAP; 256];
    let order = *b"ARNDCQEGHILKMFPSTWYV";
    let mut i = 0usize;
    while i < 20 {
        table[order[i] as usize] = i as u8;
        i += 1;
    }
    table
};

struct BitAlignment {
    n_taxa: usize,
    /// Original column count, i.e. the bootstrap sample size.
    len: usize,
    /// Retained (potentially pair-dependent) columns.
    kept: usize,
    /// u64 words per plane.
    words: usize,
    /// Taxon-major; within a taxon the six planes of a block are adjacent.
    planes: Vec<u64>,
    /// Original column -> retained index, [`CONST`] or [`DROP`].
    col_map: Vec<u32>,
    /// Number of constant gap-free columns folded away.
    const_cols: u32,
}

impl BitAlignment {
    #[inline]
    fn taxon(&self, i: usize) -> &[u64] {
        let stride = self.words * PLANES;
        &self.planes[i * stride..(i + 1) * stride]
    }

    /// Per-column validity and mismatch bitmasks for one taxon pair, computed once and
    /// then reused by every bootstrap replicate.
    #[inline]
    fn pair_masks(&self, i: usize, j: usize, valid: &mut [u64], mismatch: &mut [u64]) {
        let (a, b) = (self.taxon(i), self.taxon(j));
        for (((x, y), v), m) in a
            .chunks_exact(PLANES)
            .zip(b.chunks_exact(PLANES))
            .zip(valid.iter_mut())
            .zip(mismatch.iter_mut())
        {
            let both = x[VALID] & y[VALID];
            let diff =
                (x[0] ^ y[0]) | (x[1] ^ y[1]) | (x[2] ^ y[2]) | (x[3] ^ y[3]) | (x[4] ^ y[4]);
            *v = both;
            *m = diff & both;
        }
    }
}

fn build_bit_alignment(
    seqs: &[Vec<u8>],
    len: usize,
    pool: &ThreadPool,
) -> Result<BitAlignment, String> {
    if len > (u32::MAX - 2) as usize {
        return Err("alignment is too long for bootstrap column indexing".to_string());
    }
    let n = seqs.len();

    // Classify columns in parallel over blocks: 0 = all gap, 1 = constant and gap-free,
    // 2 = retained. `flags` bits: 1 = a gap seen, 2 = a difference seen, 4 = a residue seen.
    const BLOCK: usize = 1 << 16;
    let mut class = vec![0u8; len];
    pool.install(|| {
        class.par_chunks_mut(BLOCK).enumerate().for_each(|(b, out)| {
            let range = b * BLOCK..b * BLOCK + out.len();
            let first: Vec<u8> = seqs[0][range.clone()]
                .iter()
                .map(|&x| AA_CODE[x as usize])
                .collect();
            let mut flags: Vec<u8> = first
                .iter()
                .map(|&c| if c == GAP { 1 } else { 4 })
                .collect();
            for seq in &seqs[1..] {
                for ((&x, &r), f) in seq[range.clone()].iter().zip(&first).zip(flags.iter_mut()) {
                    let c = AA_CODE[x as usize];
                    *f |= if c == GAP { 1 } else { 4 } | (((c != r) as u8) << 1);
                }
            }
            for (o, &f) in out.iter_mut().zip(&flags) {
                *o = if f & 4 == 0 {
                    0
                } else if f & 3 == 0 {
                    1
                } else {
                    2
                };
            }
        });
    });

    let mut col_map = vec![0u32; len];
    let mut kept_cols: Vec<u32> = Vec::new();
    let mut const_cols = 0u32;
    for (column, &kind) in class.iter().enumerate() {
        col_map[column] = match kind {
            0 => DROP,
            1 => {
                const_cols += 1;
                CONST
            }
            _ => {
                kept_cols.push(column as u32);
                (kept_cols.len() - 1) as u32
            }
        };
    }

    // Bit-slice the retained columns, one task per taxon. A whole 64-column block is
    // accumulated in registers, branch-free, and stored once.
    let words = (kept_cols.len() + 63) / 64;
    let stride = words * PLANES;
    let mut planes = vec![0u64; n * stride];
    if stride > 0 {
        pool.install(|| {
            planes
                .par_chunks_mut(stride)
                .zip(seqs)
                .for_each(|(out, seq)| {
                    for (block, columns) in kept_cols.chunks(64).enumerate() {
                        let mut acc = [0u64; PLANES];
                        for (bit, &column) in columns.iter().enumerate() {
                            let code = AA_CODE[seq[column as usize] as usize];
                            let valid = (code != GAP) as u64;
                            let value = (code as u64) * valid;
                            for (plane, slot) in acc[..VALID].iter_mut().enumerate() {
                                *slot |= ((value >> plane) & 1) << bit;
                            }
                            acc[VALID] |= valid << bit;
                        }
                        out[block * PLANES..block * PLANES + PLANES].copy_from_slice(&acc);
                    }
                });
        });
    }

    Ok(BitAlignment {
        n_taxa: n,
        len,
        kept: kept_cols.len(),
        words,
        planes,
        col_map,
        const_cols,
    })
}

// ------------------------------
// Bootstrap replicate weights
// ------------------------------

/// One replicate: per-retained-column multiplicities as bit planes, plus the weight
/// that fell on constant gap-free columns.
struct Replicate {
    n_planes: usize,
    planes: Vec<u64>,
    const_weight: u32,
}

/// Draw every replicate. The RNG draws are identical to resampling `align.len` columns
/// per replicate; only the storage differs.
fn bootstrap_weights(align: &BitAlignment, replicates: usize) -> Vec<Replicate> {
    let mut rng = BootstrapRng::new(BOOTSTRAP_SEED);
    let mut weights = vec![0u32; align.kept];
    let mut out = Vec::with_capacity(replicates);

    for _ in 0..replicates {
        let mut const_weight = 0u32;
        for _ in 0..align.len {
            match align.col_map[rng.sample_index(align.len)] {
                DROP => {}
                CONST => const_weight += 1,
                slot => weights[slot as usize] += 1,
            }
        }

        let max = weights.iter().copied().max().unwrap_or(0);
        let n_planes = (u32::BITS - max.leading_zeros()) as usize;
        let mut planes = vec![0u64; n_planes * align.words];
        for (block, block_weights) in weights.chunks_mut(64).enumerate() {
            let mut acc = [0u64; u32::BITS as usize];
            for (bit, weight) in block_weights.iter_mut().enumerate() {
                for (plane, slot) in acc[..n_planes].iter_mut().enumerate() {
                    *slot |= ((*weight >> plane) as u64 & 1) << bit;
                }
                *weight = 0;
            }
            for (plane, &value) in acc[..n_planes].iter().enumerate() {
                planes[plane * align.words + block] = value;
            }
        }
        out.push(Replicate {
            n_planes,
            planes,
            const_weight,
        });
    }
    out
}

/// Weighted valid-site and mismatch counts for one pair under one replicate:
/// `Σ_c w_c·mask_c = Σ_p 2^p · popcount(mask & plane_p)`.
#[inline]
fn weighted_counts(valid: &[u64], mismatch: &[u64], replicate: &Replicate) -> (u64, u64) {
    let words = valid.len();
    let (mut valid_sites, mut mismatches) = (0u64, 0u64);
    for plane in 0..replicate.n_planes {
        let bits = &replicate.planes[plane * words..(plane + 1) * words];
        let (mut v, mut m) = (0u64, 0u64);
        for ((&w, &vm), &mm) in bits.iter().zip(valid).zip(mismatch) {
            v += (w & vm).count_ones() as u64;
            m += (w & mm).count_ones() as u64;
        }
        valid_sites += v << plane;
        mismatches += m << plane;
    }
    (valid_sites + replicate.const_weight as u64, mismatches)
}

/// Pairwise distances for one batch of replicates, as one packed row per taxon:
/// row `i` holds `(n - i - 1) * batch.len()` values, replicate-minor.
fn batch_distances(align: &BitAlignment, batch: &[Replicate], pool: &ThreadPool) -> Vec<Vec<f64>> {
    let (n, words, width) = (align.n_taxa, align.words, batch.len());
    pool.install(|| {
        (0..n)
            .into_par_iter()
            .map_init(
                || (vec![0u64; words], vec![0u64; words]),
                |(valid, mismatch), i| {
                    let mut row = vec![0.0f64; (n - i - 1) * width];
                    for j in (i + 1)..n {
                        align.pair_masks(i, j, valid, mismatch);
                        let base = (j - i - 1) * width;
                        for (slot, replicate) in row[base..base + width].iter_mut().zip(batch) {
                            let (v, m) = weighted_counts(valid, mismatch, replicate);
                            *slot = corrected_distance(v, m);
                        }
                    }
                    row
                },
            )
            .collect()
    })
}

/// The original alignment expressed as a replicate whose every column has weight one,
/// so the reference tree reuses [`batch_distances`] unchanged.
fn unit_replicate(align: &BitAlignment) -> Replicate {
    Replicate {
        n_planes: 1,
        planes: vec![u64::MAX; align.words],
        const_weight: align.const_cols,
    }
}

// ------------------------------
// Neighbor joining
// ------------------------------

/// Q-scan block width: wide enough to vectorise, narrow enough that the rare in-order
/// fix-up stays cheap.
const Q_BLOCK: usize = 8;
/// Recompute row sums exactly every this many joins, to bound incremental drift.
const ROW_SUM_REFRESH: usize = 128;

/// Order-independent key used to break exact Q ties reproducibly.
#[inline]
fn pair_key(a: u32, b: u32) -> (u32, u32) {
    (a.min(b), a.max(b))
}

/// Sum a matrix row with four accumulators to break the add dependency chain.
#[inline]
fn row_sum(row: &[f64]) -> f64 {
    let mut acc = [0.0f64; 4];
    let mut blocks = row.chunks_exact(4);
    for block in blocks.by_ref() {
        for (a, &v) in acc.iter_mut().zip(block) {
            *a += v;
        }
    }
    (acc[0] + acc[1]) + (acc[2] + acc[3]) + blocks.remainder().iter().sum::<f64>()
}

/// Per-thread scratch for one neighbor-joining run and its split bookkeeping.
///
/// `dist` is the compact n×n matrix: active nodes occupy slots `0..m`, so both O(n³)
/// scans read contiguous memory. `children[k]` are the children of internal node
/// `n + k`, the synthetic degree-2 root last.
struct Nj {
    dist: Vec<f64>,
    row_sums: Vec<f64>,
    new_row: Vec<f64>,
    slot: Vec<u32>,
    children: Vec<[u32; 2]>,
    lengths: Vec<[f64; 2]>,
    bitsets: Vec<u64>,
    split: Vec<u64>,
    stamp: Vec<u32>,
    counts: Vec<usize>,
    tag: u32,
}

impl Nj {
    fn new(n: usize, splits: usize) -> Self {
        let words = (n + 63) / 64;
        Self {
            dist: vec![0.0; n * n],
            row_sums: vec![0.0; n],
            new_row: vec![0.0; n],
            slot: Vec::with_capacity(n),
            children: Vec::with_capacity(n),
            lengths: Vec::with_capacity(n),
            bitsets: vec![0; (2 * n - 1) * words],
            split: vec![0; words],
            stamp: vec![u32::MAX; splits],
            counts: vec![0; splits],
            tag: 0,
        }
    }

    /// Copy one replicate's leaf distances out of [`batch_distances`] output.
    fn load(&mut self, n: usize, rows: &[Vec<f64>], width: usize, replicate: usize) {
        for i in 0..n {
            for (offset, j) in ((i + 1)..n).enumerate() {
                let d = rows[i][offset * width + replicate];
                self.dist[i * n + j] = d;
                self.dist[j * n + i] = d;
            }
            self.dist[i * n + i] = 0.0;
        }
    }

    /// Neighbor joining over `self.dist`, with leaf `s` initially in slot `s`.
    ///
    /// A join writes the new node into the first slot and swap-removes the second, so
    /// the active set stays packed. Row sums are maintained incrementally (O(n²)
    /// instead of O(n³)), leaving only the Q scan at O(n³).
    fn reduce(&mut self, n: usize) {
        let (dist, r, new_row, slot) = (
            &mut self.dist[..],
            &mut self.row_sums[..],
            &mut self.new_row[..],
            &mut self.slot,
        );
        self.children.clear();
        self.lengths.clear();
        slot.clear();
        slot.extend(0..n as u32);
        let mut m = n;
        for s in 0..m {
            r[s] = row_sum(&dist[s * n..s * n + m]);
        }
        let mut since_refresh = 0usize;

        while m > 2 {
            // Select the pair minimising the NJ Q criterion. The bulk of the scan is a
            // branch-free block whose minimum is compared once against the incumbent;
            // only a block that can match or beat it is re-scanned element by element.
            let factor = (m - 2) as f64;
            let mut best_q = f64::INFINITY;
            let mut best = (0usize, 1usize);
            let mut best_key = (u32::MAX, u32::MAX);
            for a in 0..m - 1 {
                let ra = r[a];
                let row = &dist[a * n..a * n + m];
                let mut b = a + 1;
                while b + Q_BLOCK <= m {
                    let mut q = [0.0f64; Q_BLOCK];
                    let mut low = f64::INFINITY;
                    for (t, value) in q.iter_mut().enumerate() {
                        *value = factor * row[b + t] - ra - r[b + t];
                        low = low.min(*value);
                    }
                    if low <= best_q {
                        for (t, &value) in q.iter().enumerate() {
                            let key = pair_key(slot[a], slot[b + t]);
                            if value < best_q || (value == best_q && key < best_key) {
                                (best_q, best, best_key) = (value, (a, b + t), key);
                            }
                        }
                    }
                    b += Q_BLOCK;
                }
                for b in b..m {
                    let value = factor * row[b] - ra - r[b];
                    let key = pair_key(slot[a], slot[b]);
                    if value < best_q || (value == best_q && key < best_key) {
                        (best_q, best, best_key) = (value, (a, b), key);
                    }
                }
            }

            let (a, b) = best;
            let dab = dist[a * n + b];
            let mut la = 0.5 * dab + (r[a] - r[b]) / (2.0 * factor);
            if la < 0.0 {
                la = 0.0;
            }
            let lb = (dab - la).max(0.0);
            self.children.push([slot[a], slot[b]]);
            self.lengths.push([la, lb]);
            slot[a] = (n + self.children.len() - 1) as u32;

            // Distances from the new node, plus the matching row-sum update
            // r[t] <- r[t] - d(a,t) - d(b,t) + d(u,t).
            let mut ra = 0.0;
            for t in 0..m {
                if t == a || t == b {
                    continue;
                }
                let (dat, dbt) = (dist[a * n + t], dist[b * n + t]);
                let dut = 0.5 * ((dat + dbt) - dab);
                new_row[t] = dut;
                r[t] = ((r[t] - dat) - dbt) + dut;
                ra += dut;
            }
            for t in 0..m {
                if t == a || t == b {
                    continue;
                }
                dist[a * n + t] = new_row[t];
                dist[t * n + a] = new_row[t];
            }
            dist[a * n + a] = 0.0;
            r[a] = ra;

            // Swap-remove slot b so the active set stays packed in 0..m-1.
            let last = m - 1;
            if b != last {
                dist.copy_within(last * n..last * n + m, b * n);
                dist[b * n + b] = 0.0;
                for t in 0..m {
                    dist[t * n + b] = dist[b * n + t];
                }
                r[b] = r[last];
                slot[b] = slot[last];
            }
            m -= 1;

            since_refresh += 1;
            if since_refresh == ROW_SUM_REFRESH {
                since_refresh = 0;
                for s in 0..m {
                    r[s] = row_sum(&dist[s * n..s * n + m]);
                }
            }
        }

        // The final two active nodes are connected under a synthetic root.
        let len = (dist[1] * 0.5).max(0.0);
        self.children.push([slot[0], slot[1]]);
        self.lengths.push([len, len]);
    }

    /// Fill descendant-leaf bitsets for every node.
    fn fill_bitsets(&mut self, n: usize) {
        let words = (n + 63) / 64;
        self.bitsets.fill(0);
        for leaf in 0..n {
            self.bitsets[leaf * words + leaf / 64] |= 1u64 << (leaf % 64);
        }
        for k in 0..self.children.len() {
            let child = self.children[k];
            let (head, tail) = self.bitsets.split_at_mut((n + k) * words);
            for (w, slot) in tail[..words].iter_mut().enumerate() {
                *slot = head[child[0] as usize * words + w] | head[child[1] as usize * words + w];
            }
        }
    }

    /// Count this tree's non-trivial splits that also occur in the reference tree.
    /// The synthetic degree-2 root makes one split appear on both of its children, so
    /// a stamp keeps each split counted once per replicate.
    fn count_splits(&mut self, n: usize, index: &HashMap<Vec<u64>, u32>) {
        let words = (n + 63) / 64;
        self.fill_bitsets(n);
        self.tag += 1;
        for node in n..n + self.children.len() {
            if !canonical_split(&self.bitsets[node * words..(node + 1) * words], n, &mut self.split)
            {
                continue;
            }
            if let Some(&id) = index.get(self.split.as_slice()) {
                let id = id as usize;
                if self.stamp[id] != self.tag {
                    self.stamp[id] = self.tag;
                    self.counts[id] += 1;
                }
            }
        }
    }
}

fn build_nodes(names: &[String], nj: &Nj) -> (Vec<Node>, usize) {
    let mut nodes: Vec<Node> = names
        .iter()
        .map(|name| Node {
            name: Some(name.clone()),
            children: Vec::new(),
        })
        .collect();
    for (child, len) in nj.children.iter().zip(&nj.lengths) {
        nodes.push(Node {
            name: None,
            children: vec![
                Edge {
                    child: child[0] as usize,
                    len: len[0],
                },
                Edge {
                    child: child[1] as usize,
                    len: len[1],
                },
            ],
        });
    }
    let root = nodes.len() - 1;
    (nodes, root)
}

/// Write the root-independent canonical form of one split into `out`; returns false for
/// trivial (terminal) splits, which carry no bootstrap information.
fn canonical_split(bits: &[u64], n_leaves: usize, out: &mut [u64]) -> bool {
    let inside: usize = bits.iter().map(|w| w.count_ones() as usize).sum();
    let outside = n_leaves - inside;
    if inside < 2 || outside < 2 {
        return false;
    }

    for (dst, &w) in out.iter_mut().zip(bits) {
        *dst = !w;
    }
    let used_last = n_leaves % 64;
    if used_last != 0 {
        let last = out.len() - 1;
        out[last] &= (1u64 << used_last) - 1;
    }
    if !(outside < inside || (outside == inside && &out[..] < bits)) {
        out.copy_from_slice(bits);
    }
    true
}

/// Index every unique non-trivial split of the reference tree.
fn index_splits(nj: &mut Nj, n: usize) -> (HashMap<Vec<u64>, u32>, Vec<Option<u32>>) {
    let words = (n + 63) / 64;
    nj.fill_bitsets(n);
    let mut index = HashMap::new();
    let mut per_node = vec![None; n + nj.children.len()];
    for node in n..n + nj.children.len() {
        if canonical_split(&nj.bitsets[node * words..(node + 1) * words], n, &mut nj.split) {
            let next = index.len() as u32;
            per_node[node] = Some(*index.entry(nj.split.clone()).or_insert(next));
        }
    }
    (index, per_node)
}

// ------------------------------
// Newick output
// ------------------------------

/// Quote/escape names containing Newick-special characters.
fn escape_name(name: &str) -> String {
    let needs_quotes = name
        .chars()
        .any(|c| matches!(c, ' ' | ':' | '(' | ')' | ',' | ';'));
    if needs_quotes {
        format!("'{}'", name.replace('\'', "''"))
    } else {
        name.to_string()
    }
}

/// Format bootstrap support as a whole-number percentage.
fn format_support(value: f64) -> String {
    format!("{}", value.round() as u32)
}

/// Append a subtree in Newick form to one shared buffer (linear in output length).
/// Internal-node labels are used for bootstrap percentages when supplied.
fn write_subtree(nodes: &[Node], node_id: usize, support: Option<&[Option<f64>]>, out: &mut String) {
    let node = &nodes[node_id];
    if node.children.is_empty() {
        out.push_str(&escape_name(node.name.as_deref().unwrap_or("")));
        return;
    }
    out.push('(');
    for (position, edge) in node.children.iter().enumerate() {
        if position > 0 {
            out.push(',');
        }
        write_subtree(nodes, edge.child, support, out);
        let _ = write!(out, ":{:.6}", edge.len);
    }
    out.push(')');
    if let Some(Some(value)) = support.and_then(|values| values.get(node_id)) {
        out.push_str(&format_support(*value));
    }
}

fn to_newick(nodes: &[Node], root: usize) -> String {
    let mut out = String::new();
    write_subtree(nodes, root, None, &mut out);
    out.push(';');
    out
}

/// Emit the bootstrap tree as an unrooted Newick tree.
///
/// `Nj::reduce` finishes by inserting a synthetic degree-2 root in the middle of the
/// last NJ edge. For an unrooted tree this root has no biological meaning, and its two
/// child branches describe complementary sides of the same split. Suppress that
/// artificial root by rooting the Newick representation at one of its internal
/// children. The output root is then a trifurcation, as expected for an unrooted binary
/// tree, and every non-trivial internal edge can carry one support label.
fn to_unrooted_newick_with_support(nodes: &[Node], root: usize, support: &[Option<f64>]) -> String {
    let mut out = String::new();
    let edges = &nodes[root].children;
    // For n >= 3 at least one side of the synthetic root is an internal node; use it as
    // the display root and reconnect the opposite side with the summed half-edges.
    let sides = if edges.len() != 2 {
        None
    } else if !nodes[edges[0].child].children.is_empty() {
        Some((edges[0].child, edges[1].child))
    } else if !nodes[edges[1].child].children.is_empty() {
        Some((edges[1].child, edges[0].child))
    } else {
        None // two-taxon tree: no internal split to support
    };

    let Some((display_root, other)) = sides else {
        write_subtree(nodes, root, Some(support), &mut out);
        out.push(';');
        return out;
    };

    out.push('(');
    for edge in &nodes[display_root].children {
        write_subtree(nodes, edge.child, Some(support), &mut out);
        let _ = write!(out, ":{:.6},", edge.len);
    }
    write_subtree(nodes, other, Some(support), &mut out);
    let _ = write!(out, ":{:.6});", edges[0].len + edges[1].len);
    out
}

// ------------------------------
// Entry points
// ------------------------------

fn validate_alignment(names: &[String], seqs: &[Vec<u8>]) -> Result<usize, String> {
    if names.len() != seqs.len() {
        return Err("names and sequences length mismatch".to_string());
    }
    if names.len() < 2 {
        return Err("need at least 2 sequences".to_string());
    }

    let alignment_len = seqs[0].len();
    if seqs.iter().any(|seq| seq.len() != alignment_len) {
        return Err("all sequences must have identical lengths".to_string());
    }
    Ok(alignment_len)
}

fn build_thread_pool(num_threads: usize) -> ThreadPool {
    ThreadPoolBuilder::new()
        .num_threads(num_threads.max(1))
        .build()
        .expect("Failed to build Rayon thread pool for the neighbor-joining stage")
}

/// Build a neighbor-joining tree and return it as Newick.
/// Requires at least two sequences of identical length.
pub fn nj_tree_newick(
    names: &[String],
    seqs: &[Vec<u8>],
    num_threads: usize,
) -> Result<String, String> {
    let len = validate_alignment(names, seqs)?;
    let n = names.len();
    let pool = build_thread_pool(num_threads);
    let align = build_bit_alignment(seqs, len, &pool)?;

    let mut nj = Nj::new(n, 0);
    nj.load(n, &batch_distances(&align, &[unit_replicate(&align)], &pool), 1, 0);
    nj.reduce(n);
    let (nodes, root) = build_nodes(names, &nj);
    Ok(to_newick(&nodes, root))
}

/// Build the NJ tree from the original alignment and annotate its internal edges with
/// non-parametric bootstrap support from column-resampled alignment replicates.
pub fn nj_tree_newick_bootstrap(
    names: &[String],
    seqs: &[Vec<u8>],
    num_threads: usize,
    replicates: usize,
) -> Result<String, String> {
    let len = validate_alignment(names, seqs)?;
    if replicates == 0 {
        return Err("bootstrap replicates must be >=1".to_string());
    }
    if len == 0 {
        return Err("cannot bootstrap an empty alignment".to_string());
    }
    let n = names.len();
    let pool = build_thread_pool(num_threads);

    // Bit-slice once; the reference tree and every replicate reuse it.
    let align = build_bit_alignment(seqs, len, &pool)?;
    let mut reference = Nj::new(n, 0);
    reference.load(n, &batch_distances(&align, &[unit_replicate(&align)], &pool), 1, 0);
    reference.reduce(n);
    let (index, per_node) = index_splits(&mut reference, n);
    let (nodes, root) = build_nodes(names, &reference);
    drop(reference);

    let weights = bootstrap_weights(&align, replicates);

    // Replicates are processed in batches: a batch's packed distances are bounded in
    // size, and its weight planes - re-read once per taxon pair - stay cache-resident.
    // The batch size never affects the result.
    let pair_count = n * (n - 1) / 2;
    let plane_words = weights.iter().map(|w| w.planes.len()).max().unwrap_or(0);
    let batch = ((8 << 20) / (plane_words * 8).max(1))
        .max(num_threads.max(1))
        .min((512 << 20) / (pair_count * 8))
        .clamp(1, replicates);

    let mut counts = vec![0usize; index.len()];
    for start in (0..replicates).step_by(batch) {
        let window = &weights[start..(start + batch).min(replicates)];
        let rows = batch_distances(&align, window, &pool);
        let width = window.len();

        // Each replicate is an independent reduction, so the whole batch runs in
        // parallel. Workspaces are recycled through a pool rather than created per
        // task: an n*n matrix is large at a few thousand taxa, and a per-task
        // workspace would spend most of the batch zeroing fresh ones.
        let workspaces: Mutex<Vec<Nj>> = Mutex::new(Vec::new());
        pool.install(|| {
            (0..width).into_par_iter().for_each(|replicate| {
                let mut nj = workspaces
                    .lock()
                    .expect("workspace pool poisoned")
                    .pop()
                    .unwrap_or_else(|| Nj::new(n, index.len()));
                nj.load(n, &rows, width, replicate);
                nj.reduce(n);
                nj.count_splits(n, &index);
                workspaces.lock().expect("workspace pool poisoned").push(nj);
            });
        });
        // Split counts are summed, so the order workspaces come back in is irrelevant.
        for nj in workspaces.into_inner().expect("workspace pool poisoned") {
            for (total, value) in counts.iter_mut().zip(nj.counts) {
                *total += value;
            }
        }
    }

    // Assign support to every node-side representation of a reference split, including
    // both complementary children of the synthetic degree-2 root. The unrooted Newick
    // writer suppresses one of those duplicates, leaving one label per internal edge.
    let mut support = vec![None; nodes.len()];
    for (node, &id) in per_node.iter().enumerate() {
        if let Some(id) = id {
            support[node] = Some(100.0 * counts[id as usize] as f64 / replicates as f64);
        }
    }
    Ok(to_unrooted_newick_with_support(&nodes, root, &support))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn labels(n: usize) -> Vec<String> {
        (0..n).map(|i| format!("t{i}")).collect()
    }

    /// Deterministic pseudo-random alignment with invariant columns, gaps and noise.
    fn synth(n: usize, len: usize, seed: u64) -> Vec<Vec<u8>> {
        let residues = b"ARNDCQEGHILKMFPSTWYV";
        let mut rng = BootstrapRng::new(seed);
        let ancestor: Vec<u8> = (0..len).map(|_| residues[rng.sample_index(20)]).collect();
        let variable: Vec<bool> = (0..len).map(|_| rng.sample_index(100) < 60).collect();
        (0..n)
            .map(|_| {
                (0..len)
                    .map(|c| match rng.sample_index(100) {
                        _ if !variable[c] => ancestor[c],
                        d if d < 5 => b'-',
                        d if d < 8 => b'X',
                        d if d < 40 => residues[rng.sample_index(20)],
                        _ => ancestor[c],
                    })
                    .collect()
            })
            .collect()
    }

    fn naive_distance(a: &[u8], b: &[u8], weights: Option<&[u32]>) -> f64 {
        let (mut valid, mut mismatches) = (0u64, 0u64);
        for (column, (&x, &y)) in a.iter().zip(b).enumerate() {
            let (cx, cy) = (AA_CODE[x as usize], AA_CODE[y as usize]);
            let w = weights.map_or(1, |w| w[column] as u64);
            if cx < GAP && cy < GAP {
                valid += w;
                if cx != cy {
                    mismatches += w;
                }
            }
        }
        corrected_distance(valid, mismatches)
    }

    /// Textbook NJ over a (2n-1)² matrix with an explicit active list, i.e. the
    /// previous implementation's structure, returning splits only.
    fn reference_nj(n: usize, leaf: &[f64]) -> std::collections::BTreeSet<Vec<u64>> {
        let dim = 2 * n - 1;
        let words = (n + 63) / 64;
        let mut dist = vec![0.0f64; dim * dim];
        let mut bits = vec![0u64; dim * words];
        for i in 0..n {
            bits[i * words + i / 64] |= 1u64 << (i % 64);
            for j in 0..n {
                dist[i * dim + j] = leaf[i * n + j];
            }
        }
        let mut active: Vec<usize> = (0..n).collect();
        let mut r = vec![0.0f64; dim];
        let mut next = n;
        while active.len() > 2 {
            let m = active.len();
            for &i in &active {
                r[i] = active.iter().map(|&j| dist[i * dim + j]).sum();
            }
            let mut best = (active[0], active[1]);
            let mut best_q = f64::INFINITY;
            let mut best_key = (usize::MAX, usize::MAX);
            for (pos, &i) in active.iter().enumerate() {
                for &j in active.iter().skip(pos + 1) {
                    let q = ((m - 2) as f64) * dist[i * dim + j] - r[i] - r[j];
                    let key = (i.min(j), i.max(j));
                    if q < best_q || (q == best_q && key < best_key) {
                        (best_q, best, best_key) = (q, (i, j), key);
                    }
                }
            }
            let (i, j) = best;
            let dij = dist[i * dim + j];
            let u = next;
            next += 1;
            for w in 0..words {
                bits[u * words + w] = bits[i * words + w] | bits[j * words + w];
            }
            for &k in &active {
                if k == i || k == j {
                    continue;
                }
                let duk = 0.5 * (dist[i * dim + k] + dist[j * dim + k] - dij);
                dist[u * dim + k] = duk;
                dist[k * dim + u] = duk;
            }
            active.retain(|&x| x != i && x != j);
            active.push(u);
        }
        let mut out = std::collections::BTreeSet::new();
        let mut split = vec![0u64; words];
        for node in n..next {
            if canonical_split(&bits[node * words..(node + 1) * words], n, &mut split) {
                out.insert(split.clone());
            }
        }
        out
    }

    fn splits_of(nj: &mut Nj, n: usize) -> std::collections::BTreeSet<Vec<u64>> {
        let words = (n + 63) / 64;
        nj.fill_bitsets(n);
        let mut out = std::collections::BTreeSet::new();
        let mut split = vec![0u64; words];
        for node in n..n + nj.children.len() {
            if canonical_split(&nj.bitsets[node * words..(node + 1) * words], n, &mut split) {
                out.insert(split.clone());
            }
        }
        out
    }

    #[test]
    fn residue_table_matches_lg_order() {
        for (index, &byte) in b"ARNDCQEGHILKMFPSTWYV".iter().enumerate() {
            assert_eq!(AA_CODE[byte as usize], index as u8);
        }
        for &byte in b"-X*?BZJUOax" {
            assert_eq!(AA_CODE[byte as usize], GAP);
        }
    }

    #[test]
    fn column_reduction_preserves_distances() {
        let seqs = synth(17, 501, 11);
        let pool = build_thread_pool(2);
        let align = build_bit_alignment(&seqs, 501, &pool).unwrap();
        assert!(align.kept < 501, "expected columns to be folded away");
        assert_eq!(
            align.kept + align.const_cols as usize
                + align.col_map.iter().filter(|&&s| s == DROP).count(),
            501
        );

        let rows = batch_distances(&align, &[unit_replicate(&align)], &pool);
        for i in 0..17 {
            for j in (i + 1)..17 {
                let expected = naive_distance(&seqs[i], &seqs[j], None);
                assert!((rows[i][j - i - 1] - expected).abs() < 1e-12, "pair ({i},{j})");
            }
        }
    }

    #[test]
    fn bootstrap_distances_match_dense_reference() {
        // Gaps, unknowns, invariant columns and an all-gap column.
        let seqs = vec![
            b"AARN-XVA-KKMAAV".to_vec(),
            b"ARRNAXVA-KKLAAV".to_vec(),
            b"NNRNA-VR-KKMAAV".to_vec(),
            b"NNRDA-VR-KKMAAC".to_vec(),
        ];
        let (n, len) = (seqs.len(), seqs[0].len());
        let pool = build_thread_pool(2);
        let align = build_bit_alignment(&seqs, len, &pool).unwrap();
        assert!(align.const_cols > 0 && align.kept < len);

        let replicates = bootstrap_weights(&align, 9);
        let rows = batch_distances(&align, &replicates, &pool);

        // Replay the same draws over the *original* columns, bypassing the reduction.
        let mut rng = BootstrapRng::new(BOOTSTRAP_SEED);
        for replicate in 0..replicates.len() {
            let mut weights = vec![0u32; len];
            for _ in 0..len {
                weights[rng.sample_index(len)] += 1;
            }
            for i in 0..n {
                for j in (i + 1)..n {
                    let expected = naive_distance(&seqs[i], &seqs[j], Some(&weights));
                    let observed = rows[i][(j - i - 1) * replicates.len() + replicate];
                    assert!(
                        (expected - observed).abs() < 1e-12,
                        "replicate {replicate} pair ({i},{j}): {observed} vs {expected}"
                    );
                }
            }
        }
    }

    #[test]
    fn compact_nj_matches_textbook_nj() {
        // Continuous random distances: exact Q ties have probability zero, so the
        // compact reduction must agree with the textbook one exactly.
        for &n in &[5usize, 8, 23, 64, 129, 200] {
            let mut rng = BootstrapRng::new(0xC0FFEE + n as u64);
            let mut dist = vec![0.0f64; n * n];
            for i in 0..n {
                for j in (i + 1)..n {
                    let d = 0.01 + (rng.next_u64() >> 11) as f64 / (1u64 << 53) as f64;
                    dist[i * n + j] = d;
                    dist[j * n + i] = d;
                }
            }
            let mut nj = Nj::new(n, 0);
            nj.dist.copy_from_slice(&dist);
            nj.reduce(n);
            assert_eq!(splits_of(&mut nj, n), reference_nj(n, &dist), "n={n}");
        }
    }

    #[test]
    fn complementary_sides_have_the_same_split_key() {
        let (mut a, mut b) = (vec![0u64; 1], vec![0u64; 1]);
        assert!(canonical_split(&[0b0011u64], 4, &mut a));
        assert!(canonical_split(&[0b1100u64], 4, &mut b));
        assert_eq!(a, b);
        assert!(!canonical_split(&[0b0001u64], 4, &mut a));
        assert!(!canonical_split(&[0b0111u64], 4, &mut a));
    }

    #[test]
    fn bootstrap_support_is_rounded_to_integer_percent() {
        assert_eq!(format_support(89.5), "90");
        assert_eq!(format_support(89.4), "89");
        assert_eq!(format_support(100.0), "100");
    }

    #[test]
    fn bootstrap_support_is_written_and_reproducible() {
        let names = vec!["A", "B", "C", "D"]
            .into_iter()
            .map(str::to_string)
            .collect::<Vec<_>>();
        let seqs = vec![
            b"AAAAAAAA".to_vec(),
            b"AAAAAAAA".to_vec(),
            b"RRRRRRRR".to_vec(),
            b"RRRRRRRR".to_vec(),
        ];
        let tree1 = nj_tree_newick_bootstrap(&names, &seqs, 1, 20).unwrap();
        let tree2 = nj_tree_newick_bootstrap(&names, &seqs, 1, 20).unwrap();
        assert_eq!(tree1, tree2);
        assert!(tree1.contains(")100:"), "unexpected bootstrap tree: {tree1}");
        assert_eq!(tree1.matches(")100:").count(), 1, "tree: {tree1}");
    }

    #[test]
    fn output_is_independent_of_thread_count() {
        let seqs = synth(24, 900, 42);
        let names = labels(24);
        let one = nj_tree_newick_bootstrap(&names, &seqs, 1, 30).unwrap();
        for threads in [2usize, 3, 7] {
            assert_eq!(
                one,
                nj_tree_newick_bootstrap(&names, &seqs, threads, 30).unwrap()
            );
        }
        assert_eq!(
            nj_tree_newick(&names, &seqs, 1).unwrap(),
            nj_tree_newick(&names, &seqs, 4).unwrap()
        );
    }
}
