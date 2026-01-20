use std::collections::HashSet;

use rayon::prelude::*;
use rayon::ThreadPoolBuilder;

#[derive(Debug, Clone)]
pub struct Edge {
    pub child: usize,
    pub len: f64,
}

#[derive(Debug, Clone)]
pub struct Node {
    pub name: Option<String>,
    pub children: Vec<Edge>,
}

// LG stationary frequencies from Le & Gascuel (2008).
const LG_PI: [f64; 20] = [
    0.079_066, 0.055_941, 0.041_977, 0.053_052, 0.012_937, 0.040_767, 0.071_586,
    0.057_337, 0.022_355, 0.062_157, 0.099_081, 0.064_600, 0.022_951, 0.042_302,
    0.044_040, 0.061_197, 0.053_287, 0.012_777, 0.027_843, 0.070_200,
];

/// Map an uppercase amino-acid byte to its LG index (A,R,N,...,V).
/// Returns None for gaps/unknowns so they are ignored in distance estimates.
fn aa_index(b: u8) -> Option<usize> {
    match b {
        b'A' => Some(0),
        b'R' => Some(1),
        b'N' => Some(2),
        b'D' => Some(3),
        b'C' => Some(4),
        b'Q' => Some(5),
        b'E' => Some(6),
        b'G' => Some(7),
        b'H' => Some(8),
        b'I' => Some(9),
        b'L' => Some(10),
        b'K' => Some(11),
        b'M' => Some(12),
        b'F' => Some(13),
        b'P' => Some(14),
        b'S' => Some(15),
        b'T' => Some(16),
        b'W' => Some(17),
        b'Y' => Some(18),
        b'V' => Some(19),
        _ => None,
    }
}

/// Compute an LG+F81 corrected distance between two aligned sequences.
/// Assumes equal-length, uppercase amino-acid bytes; ignores unknowns.
fn lg_f81_distance(a: &[u8], b: &[u8]) -> f64 {
    let mut valid_sites = 0usize;
    let mut mismatches = 0usize;

    for (&aa, &bb) in a.iter().zip(b.iter()) {
        let ia = aa_index(aa);
        let ib = aa_index(bb);
        if let (Some(ia), Some(ib)) = (ia, ib) {
            valid_sites += 1;
            if ia != ib {
                mismatches += 1;
            }
        }
    }

    if valid_sites == 0 {
        return 0.0;
    }

    let p = (mismatches as f64) / (valid_sites as f64);
    let c: f64 = 1.0 - LG_PI.iter().map(|pi| pi * pi).sum::<f64>();
    let mut pc = p / c;
    if pc >= 1.0 {
        pc = 0.999_999_999;
    }
    -(1.0 - pc).ln()
}

fn idx(i: usize, j: usize, dim: usize) -> usize {
    i * dim + j
}

/// Build the full pairwise distance matrix (size (2n-1)^2) for NJ.
fn compute_distance_matrix(seqs: &[Vec<u8>], num_threads: usize) -> Vec<f64> {
    let n = seqs.len();
    let dim = 2 * n - 1;
    let mut dist = vec![0.0; dim * dim];
    let dist_ptr = dist.as_mut_ptr() as usize;
    let n_threads = num_threads.max(1);
    let pool = ThreadPoolBuilder::new()
        .num_threads(n_threads)
        .build()
        .expect("Failed to build Rayon thread pool in compute_distance_matrix");
    pool.install(|| {
        (0..n).into_par_iter().for_each(|i| {
            for j in (i + 1)..n {
                let d = lg_f81_distance(&seqs[i], &seqs[j]);
                let ij = idx(i, j, dim);
                let ji = idx(j, i, dim);
                unsafe {
                    // Safety: each (i, j) pair is unique per thread, so writes do not overlap.
                    let dist_ptr = dist_ptr as *mut f64;
                    *dist_ptr.add(ij) = d;
                    *dist_ptr.add(ji) = d;
                }
            }
        });
    });

    dist
}

/// Perform neighbor-joining using the precomputed distances.
/// Returns the full node list and the root index.
fn neighbor_joining(names: &[String], dist: &mut [f64]) -> (Vec<Node>, usize) {
    let n = names.len();
    let dim = 2 * n - 1;
    let mut nodes = Vec::with_capacity(dim);
    for name in names {
        nodes.push(Node {
            name: Some(name.clone()),
            children: Vec::new(),
        });
    }

    let mut active: Vec<usize> = (0..n).collect();
    while active.len() > 2 {
        let m = active.len();
        let mut r = vec![0.0; dim];
        for &i in &active {
            let mut sum = 0.0;
            for &j in &active {
                if i != j {
                    sum += dist[idx(i, j, dim)];
                }
            }
            r[i] = sum;
        }

        let mut best_pair = (active[0], active[1]);
        let mut best_q = f64::INFINITY;
        for (pos_i, &i) in active.iter().enumerate() {
            for &j in active.iter().skip(pos_i + 1) {
                let q = ((m - 2) as f64) * dist[idx(i, j, dim)] - r[i] - r[j];
                if q < best_q {
                    best_q = q;
                    best_pair = (i, j);
                }
            }
        }

        let (i, j) = best_pair;
        let dij = dist[idx(i, j, dim)];
        let denom = 2.0 * ((m - 2) as f64);
        let mut li = 0.5 * dij + (r[i] - r[j]) / denom;
        let mut lj = dij - li;
        if li < 0.0 {
            li = 0.0;
        }
        if lj < 0.0 {
            lj = 0.0;
        }

        let u = nodes.len();
        nodes.push(Node {
            name: None,
            children: vec![Edge { child: i, len: li }, Edge { child: j, len: lj }],
        });

        for &k in &active {
            if k == i || k == j {
                continue;
            }
            let duk = 0.5 * (dist[idx(i, k, dim)] + dist[idx(j, k, dim)] - dij);
            dist[idx(u, k, dim)] = duk;
            dist[idx(k, u, dim)] = duk;
        }

        active.retain(|&x| x != i && x != j);
        active.push(u);
    }

    let a = active[0];
    let b = active[1];
    let mut len = dist[idx(a, b, dim)] * 0.5;
    if len < 0.0 {
        len = 0.0;
    }
    let root = nodes.len();
    nodes.push(Node {
        name: None,
        children: vec![Edge { child: a, len }, Edge { child: b, len }],
    });

    (nodes, root)
}

/// Quote/escape names containing Newick-special characters.
fn escape_name(name: &str) -> String {
    let needs_quotes = name
        .chars()
        .any(|c| matches!(c, ' ' | ':' | '(' | ')' | ',' | ';'));
    if needs_quotes {
        let escaped = name.replace('\'', "''");
        format!("'{}'", escaped)
    } else {
        name.to_string()
    }
}

/// Format branch lengths with fixed precision for stable output.
fn format_len(len: f64) -> String {
    format!("{:.6}", len)
}

/// Recursively format a subtree in Newick (without trailing semicolon).
fn format_subtree(nodes: &[Node], node_id: usize) -> String {
    let node = &nodes[node_id];
    if node.children.is_empty() {
        return escape_name(node.name.as_deref().unwrap_or(""));
    }

    let mut parts = Vec::with_capacity(node.children.len());
    for edge in &node.children {
        let child_str = format_subtree(nodes, edge.child);
        parts.push(format!("{}:{}", child_str, format_len(edge.len)));
    }
    format!("({})", parts.join(","))
}

/// Emit a full Newick string from the rooted node list.
fn to_newick(nodes: &[Node], root: usize) -> String {
    format!("{};", format_subtree(nodes, root))
}

/// Build a neighbor-joining tree and return it as Newick.
/// Requires at least two sequences of identical length.
pub fn nj_tree_newick(
    names: &[String],
    seqs: &[Vec<u8>],
    num_threads: usize,
) -> Result<String, String> {
    if names.len() != seqs.len() {
        return Err("names and sequences length mismatch".to_string());
    }
    if names.len() < 2 {
        return Err("need at least 2 sequences".to_string());
    }
    let mut lengths = HashSet::new();
    for seq in seqs {
        lengths.insert(seq.len());
    }
    if lengths.len() != 1 {
        return Err("all sequences must have identical lengths".to_string());
    }

    let mut dist = compute_distance_matrix(seqs, num_threads);
    let (nodes, root) = neighbor_joining(names, &mut dist);
    Ok(to_newick(&nodes, root))
}
