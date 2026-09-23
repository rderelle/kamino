//! Minimal neighbor-joining implementation used for the optional `--nj` tree.
//!
//! Distances are computed from pairwise amino-acid differences while ignoring gaps
//! and unknown residues. The alignment is bit-sliced so each pair comparison handles
//! 64 columns at a time, and invariant gap-free/all-gap columns are folded away.
//! The NJ reduction uses a compact matrix and incrementally maintained row sums.
use rayon::prelude::*;
use rayon::{ThreadPool, ThreadPoolBuilder};
use std::fmt::Write as _;

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

const LG_PI: [f64; 20] = [
    0.079_066, 0.055_941, 0.041_977, 0.053_052, 0.012_937, 0.040_767, 0.071_586, 0.057_337,
    0.022_355, 0.062_157, 0.099_081, 0.064_600, 0.022_951, 0.042_302, 0.044_040, 0.061_197,
    0.053_287, 0.012_777, 0.027_843, 0.070_200,
];

const C: f64 = {
    let mut sum = 0.0;
    let mut i = 0;
    while i < LG_PI.len() {
        sum += LG_PI[i] * LG_PI[i];
        i += 1;
    }
    1.0 - sum
};

const GAP: u8 = 20;
const PLANES: usize = 6;
const VALID: usize = 5;

const AA_CODE: [u8; 256] = {
    let mut table = [GAP; 256];
    let order = *b"ARNDCQEGHILKMFPSTWYV";
    let mut i = 0;
    while i < order.len() {
        table[order[i] as usize] = i as u8;
        i += 1;
    }
    table
};

fn corrected_distance(valid_sites: u64, mismatches: u64) -> f64 {
    if valid_sites == 0 {
        return 0.0;
    }
    let p = mismatches as f64 / valid_sites as f64;
    let pc = (p / C).min(0.999_999_999);
    -(1.0 - pc).ln()
}

struct BitAlignment {
    n_taxa: usize,
    words: usize,
    planes: Vec<u64>,
    const_cols: u64,
}

impl BitAlignment {
    #[inline]
    fn taxon(&self, i: usize) -> &[u64] {
        let stride = self.words * PLANES;
        &self.planes[i * stride..(i + 1) * stride]
    }

    #[inline]
    fn counts(&self, i: usize, j: usize) -> (u64, u64) {
        let (a, b) = (self.taxon(i), self.taxon(j));
        let (mut valid_sites, mut mismatches) = (self.const_cols, 0u64);
        for (x, y) in a.chunks_exact(PLANES).zip(b.chunks_exact(PLANES)) {
            let both = x[VALID] & y[VALID];
            let diff =
                (x[0] ^ y[0]) | (x[1] ^ y[1]) | (x[2] ^ y[2]) | (x[3] ^ y[3]) | (x[4] ^ y[4]);
            valid_sites += both.count_ones() as u64;
            mismatches += (diff & both).count_ones() as u64;
        }
        (valid_sites, mismatches)
    }
}

fn build_bit_alignment(seqs: &[Vec<u8>], len: usize, pool: &ThreadPool) -> BitAlignment {
    // 0 = all gap, 1 = constant and gap-free, 2 = potentially pair-dependent.
    const BLOCK: usize = 1 << 16;
    let mut class = vec![0u8; len];
    pool.install(|| {
        class
            .par_chunks_mut(BLOCK)
            .enumerate()
            .for_each(|(block, out)| {
                let range = block * BLOCK..block * BLOCK + out.len();
                let first: Vec<u8> = seqs[0][range.clone()]
                    .iter()
                    .map(|&byte| AA_CODE[byte as usize])
                    .collect();
                // Flags: 1 = gap, 2 = difference, 4 = residue.
                let mut flags: Vec<u8> = first
                    .iter()
                    .map(|&code| if code == GAP { 1 } else { 4 })
                    .collect();
                for seq in &seqs[1..] {
                    for ((&byte, &reference), flag) in
                        seq[range.clone()].iter().zip(&first).zip(&mut flags)
                    {
                        let code = AA_CODE[byte as usize];
                        *flag |= if code == GAP {
                            1
                        } else {
                            4 | (((code != reference) as u8) << 1)
                        };
                    }
                }
                for (kind, &flag) in out.iter_mut().zip(&flags) {
                    *kind = if flag & 4 == 0 {
                        0
                    } else if flag & 3 == 0 {
                        1
                    } else {
                        2
                    };
                }
            });
    });

    let kept_cols: Vec<usize> = class
        .iter()
        .enumerate()
        .filter_map(|(column, &kind)| (kind == 2).then_some(column))
        .collect();
    let const_cols = class.iter().filter(|&&kind| kind == 1).count() as u64;
    let words = kept_cols.len().div_ceil(64);
    let stride = words * PLANES;
    let mut planes = vec![0u64; seqs.len() * stride];
    if stride != 0 {
        pool.install(|| {
            planes
                .par_chunks_mut(stride)
                .zip(seqs)
                .for_each(|(out, seq)| {
                    for (block, columns) in kept_cols.chunks(64).enumerate() {
                        let mut acc = [0u64; PLANES];
                        for (bit, &column) in columns.iter().enumerate() {
                            let code = AA_CODE[seq[column] as usize];
                            let valid = (code != GAP) as u64;
                            let value = code as u64 * valid;
                            for (plane, slot) in acc[..VALID].iter_mut().enumerate() {
                                *slot |= ((value >> plane) & 1) << bit;
                            }
                            acc[VALID] |= valid << bit;
                        }
                        out[block * PLANES..(block + 1) * PLANES].copy_from_slice(&acc);
                    }
                });
        });
    }
    BitAlignment {
        n_taxa: seqs.len(),
        words,
        planes,
        const_cols,
    }
}

fn compute_distance_matrix(align: &BitAlignment, pool: &ThreadPool) -> Vec<Vec<f64>> {
    pool.install(|| {
        (0..align.n_taxa)
            .into_par_iter()
            .map(|i| {
                ((i + 1)..align.n_taxa)
                    .map(|j| {
                        let (valid, mismatches) = align.counts(i, j);
                        corrected_distance(valid, mismatches)
                    })
                    .collect()
            })
            .collect()
    })
}

const Q_BLOCK: usize = 8;
const ROW_SUM_REFRESH: usize = 128;

#[inline]
fn row_sum(row: &[f64]) -> f64 {
    let mut acc = [0.0; 4];
    let mut blocks = row.chunks_exact(4);
    for block in blocks.by_ref() {
        for (sum, &value) in acc.iter_mut().zip(block) {
            *sum += value;
        }
    }
    (acc[0] + acc[1]) + (acc[2] + acc[3]) + blocks.remainder().iter().sum::<f64>()
}

struct Nj {
    dist: Vec<f64>,
    row_sums: Vec<f64>,
    new_row: Vec<f64>,
    slot: Vec<u32>,
    children: Vec<[u32; 2]>,
    lengths: Vec<[f64; 2]>,
}

impl Nj {
    fn new(n: usize) -> Self {
        Self {
            dist: vec![0.0; n * n],
            row_sums: vec![0.0; n],
            new_row: vec![0.0; n],
            slot: Vec::with_capacity(n),
            children: Vec::with_capacity(n),
            lengths: Vec::with_capacity(n),
        }
    }

    fn load(&mut self, n: usize, rows: &[Vec<f64>]) {
        for (i, row) in rows.iter().enumerate().take(n) {
            for (offset, j) in ((i + 1)..n).enumerate() {
                let distance = row[offset];
                self.dist[i * n + j] = distance;
                self.dist[j * n + i] = distance;
            }
        }
    }

    fn reduce(&mut self, n: usize) {
        let (dist, sums, new_row, slot) = (
            &mut self.dist[..],
            &mut self.row_sums[..],
            &mut self.new_row[..],
            &mut self.slot,
        );
        slot.extend(0..n as u32);
        let mut m = n;
        for i in 0..m {
            sums[i] = row_sum(&dist[i * n..i * n + m]);
        }
        let mut since_refresh = 0;

        while m > 2 {
            let factor = (m - 2) as f64;
            let mut best_q = f64::INFINITY;
            let mut best = (0, 1);
            let mut best_key = (u32::MAX, u32::MAX);
            for a in 0..m - 1 {
                let row = &dist[a * n..a * n + m];
                let mut b = a + 1;
                while b + Q_BLOCK <= m {
                    let mut q = [0.0; Q_BLOCK];
                    let mut low = f64::INFINITY;
                    for (offset, value) in q.iter_mut().enumerate() {
                        *value = factor * row[b + offset] - sums[a] - sums[b + offset];
                        low = low.min(*value);
                    }
                    if low <= best_q {
                        for (offset, &value) in q.iter().enumerate() {
                            let key = pair_key(slot[a], slot[b + offset]);
                            if value < best_q || (value == best_q && key < best_key) {
                                (best_q, best, best_key) = (value, (a, b + offset), key);
                            }
                        }
                    }
                    b += Q_BLOCK;
                }
                for b in b..m {
                    let value = factor * row[b] - sums[a] - sums[b];
                    let key = pair_key(slot[a], slot[b]);
                    if value < best_q || (value == best_q && key < best_key) {
                        (best_q, best, best_key) = (value, (a, b), key);
                    }
                }
            }

            let (a, b) = best;
            let dab = dist[a * n + b];
            let la = (0.5 * dab + (sums[a] - sums[b]) / (2.0 * factor)).max(0.0);
            let lb = (dab - la).max(0.0);
            self.children.push([slot[a], slot[b]]);
            self.lengths.push([la, lb]);
            slot[a] = (n + self.children.len() - 1) as u32;

            let mut new_sum = 0.0;
            for t in 0..m {
                if t != a && t != b {
                    let (dat, dbt) = (dist[a * n + t], dist[b * n + t]);
                    let dut = 0.5 * (dat + dbt - dab);
                    new_row[t] = dut;
                    sums[t] = sums[t] - dat - dbt + dut;
                    new_sum += dut;
                }
            }
            for t in 0..m {
                if t != a && t != b {
                    dist[a * n + t] = new_row[t];
                    dist[t * n + a] = new_row[t];
                }
            }
            dist[a * n + a] = 0.0;
            sums[a] = new_sum;

            let last = m - 1;
            if b != last {
                dist.copy_within(last * n..last * n + m, b * n);
                dist[b * n + b] = 0.0;
                for t in 0..m {
                    dist[t * n + b] = dist[b * n + t];
                }
                sums[b] = sums[last];
                slot[b] = slot[last];
            }
            m -= 1;
            since_refresh += 1;
            if since_refresh == ROW_SUM_REFRESH {
                since_refresh = 0;
                for i in 0..m {
                    sums[i] = row_sum(&dist[i * n..i * n + m]);
                }
            }
        }

        let len = (dist[1] * 0.5).max(0.0);
        self.children.push([slot[0], slot[1]]);
        self.lengths.push([len, len]);
    }
}

#[inline]
fn pair_key(a: u32, b: u32) -> (u32, u32) {
    (a.min(b), a.max(b))
}

fn build_nodes(names: &[String], nj: &Nj) -> (Vec<Node>, usize) {
    let mut nodes: Vec<Node> = names
        .iter()
        .map(|name| Node {
            name: Some(name.clone()),
            children: Vec::new(),
        })
        .collect();
    for (children, lengths) in nj.children.iter().zip(&nj.lengths) {
        nodes.push(Node {
            name: None,
            children: vec![
                Edge {
                    child: children[0] as usize,
                    len: lengths[0],
                },
                Edge {
                    child: children[1] as usize,
                    len: lengths[1],
                },
            ],
        });
    }
    let root = nodes.len() - 1;
    (nodes, root)
}

fn escape_name(name: &str) -> String {
    if name
        .chars()
        .any(|character| matches!(character, ' ' | ':' | '(' | ')' | ',' | ';'))
    {
        format!("'{}'", name.replace('\'', "''"))
    } else {
        name.to_string()
    }
}

fn write_subtree(nodes: &[Node], node_id: usize, out: &mut String) {
    let node = &nodes[node_id];
    if node.children.is_empty() {
        out.push_str(&escape_name(node.name.as_deref().unwrap_or("")));
        return;
    }
    out.push('(');
    for (position, edge) in node.children.iter().enumerate() {
        if position != 0 {
            out.push(',');
        }
        write_subtree(nodes, edge.child, out);
        let _ = write!(out, ":{:.6}", edge.len);
    }
    out.push(')');
}

fn to_newick(nodes: &[Node], root: usize) -> String {
    let mut out = String::new();
    write_subtree(nodes, root, &mut out);
    out.push(';');
    out
}

fn validate_alignment(names: &[String], seqs: &[Vec<u8>]) -> Result<usize, String> {
    if names.len() != seqs.len() {
        return Err("names and sequences length mismatch".to_string());
    }
    if names.len() < 2 {
        return Err("need at least 2 sequences".to_string());
    }
    let len = seqs[0].len();
    if seqs.iter().any(|seq| seq.len() != len) {
        return Err("all sequences must have identical lengths".to_string());
    }
    Ok(len)
}

pub fn nj_tree_newick(
    names: &[String],
    seqs: &[Vec<u8>],
    num_threads: usize,
) -> Result<String, String> {
    let len = validate_alignment(names, seqs)?;
    let pool = ThreadPoolBuilder::new()
        .num_threads(num_threads.max(1))
        .build()
        .expect("Failed to build Rayon thread pool for the neighbor-joining stage");
    let alignment = build_bit_alignment(seqs, len, &pool);
    let mut nj = Nj::new(names.len());
    nj.load(names.len(), &compute_distance_matrix(&alignment, &pool));
    nj.reduce(names.len());
    let (nodes, root) = build_nodes(names, &nj);
    Ok(to_newick(&nodes, root))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn naive_distance(a: &[u8], b: &[u8]) -> f64 {
        let (mut valid, mut mismatches) = (0, 0);
        for (&a, &b) in a.iter().zip(b) {
            let (a, b) = (AA_CODE[a as usize], AA_CODE[b as usize]);
            if a < GAP && b < GAP {
                valid += 1;
                mismatches += (a != b) as u64;
            }
        }
        corrected_distance(valid, mismatches)
    }

    #[test]
    fn bit_slicing_preserves_pairwise_distances() {
        let seqs = vec![
            b"AARN-XVA-KKMAAV".to_vec(),
            b"ARRNAXVA-KKLAAV".to_vec(),
            b"NNRNA-VR-KKMAAV".to_vec(),
            b"NNRDA-VR-KKMAAC".to_vec(),
        ];
        let pool = ThreadPoolBuilder::new().num_threads(2).build().unwrap();
        let alignment = build_bit_alignment(&seqs, seqs[0].len(), &pool);
        let rows = compute_distance_matrix(&alignment, &pool);
        for i in 0..seqs.len() {
            for j in i + 1..seqs.len() {
                assert!((rows[i][j - i - 1] - naive_distance(&seqs[i], &seqs[j])).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn output_is_independent_of_thread_count() {
        let names = ["A", "B", "C", "D"].map(str::to_string);
        let seqs = vec![
            b"AAAA----ARND".to_vec(),
            b"AAAR----ARND".to_vec(),
            b"RRRR----ARNE".to_vec(),
            b"RRRN----ARNE".to_vec(),
        ];
        assert_eq!(
            nj_tree_newick(&names, &seqs, 1).unwrap(),
            nj_tree_newick(&names, &seqs, 4).unwrap()
        );
    }
}
