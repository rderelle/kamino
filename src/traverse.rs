use crate::graph::Graph;
use crate::recode::RECODE_BITS_PER_SYMBOL;
use ahash::RandomState;
use bitvec::prelude::*;
use hashbrown::{HashMap, HashSet};
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;

/// One path with the species set accumulated through bifurcations
#[derive(Clone)]
pub struct PathRec {
    pub nodes: Vec<u64>,
    pub species: BitVec<u32, Lsb0>,
}

/// VariantGroups: map (start,end) -> list of paths (with species sets)
pub type VariantGroups = HashMap<(u64, u64), Vec<PathRec>, RandomState>;

type ScoredGroup = (((u64, u64), Vec<PathRec>), f64, usize);

/// Public entry point expected by main.rs:
/// find_variant_groups(&g, &start_kmers, &end_kmers, max_depth, bubble_ratio, num_threads)
pub fn find_variant_groups(
    g: &Graph,
    start_kmers: &HashSet<u64>,
    end_kmers: &HashSet<u64>,
    max_depth: usize,
    bubble_ratio: f32,
    num_threads: usize,
) -> VariantGroups {
    // Infer number of species from node bitvectors
    let n_species = g.n_species;

    // Minimum number of species required
    let need = ((bubble_ratio as f64) * (n_species as f64)).ceil() as usize;

    // Use mask from the graph
    let mask_k1 = g.k1_mask;

    // Minimum path length (edges) before storing a variant group entry.
    let min_path_edges = g.k;

    collect_variant_groups(
        g,
        start_kmers,
        end_kmers,
        need,
        max_depth,
        mask_k1,
        min_path_edges,
        num_threads,
        n_species,
    )
}

#[allow(clippy::too_many_arguments)]
fn collect_variant_groups(
    g: &Graph,
    start_kmers: &HashSet<u64>,
    end_kmers: &HashSet<u64>,
    need: usize,
    max_depth: usize,
    mask_k1: u64,
    min_path_edges: usize,
    num_threads: usize,
    n_species: usize,
) -> VariantGroups {
    // Deterministic order for reproducibility
    let mut starts_sorted: Vec<u64> = start_kmers.iter().copied().collect();
    starts_sorted.sort_unstable();

    // Build a local Rayon thread pool with the requested number of threads
    let n_threads = num_threads.max(1);
    let pool = ThreadPoolBuilder::new()
        .num_threads(n_threads)
        .build()
        .expect("Failed to build Rayon thread pool in find_variant_groups");

    pool.install(|| {
        // Parallel over start nodes
        let partial_groups: Vec<VariantGroups> = starts_sorted
            .par_iter()
            .map(|&start| {
                // Per-start traversal
                let mut per_start: VariantGroups = VariantGroups::default();
                collect_paths_iterative_species_from_last_bifurcation(
                    g,
                    start,
                    end_kmers,
                    max_depth,
                    mask_k1,
                    min_path_edges,
                    &mut per_start,
                );

                // per (start,end) processing and filtering
                let mut keys: Vec<(u64, u64)> = per_start.keys().copied().collect();
                keys.sort_unstable();

                let mut kept_local: VariantGroups = VariantGroups::default();
                for key in keys {
                    if let Some(paths) = per_start.remove(&key) {
                        if let Some((k, v)) =
                            filter_one_group_3steps(g, key, paths, need, n_species)
                        {
                            kept_local.insert(k, v);
                        }
                    }
                }

                kept_local
            })
            .collect();

        // Merge all per-start results
        let mut all_kept: VariantGroups = VariantGroups::default();
        for vg in partial_groups {
            for (k, v) in vg {
                all_kept.insert(k, v);
            }
        }

        // Final global selection
        post_select_variant_groups(all_kept)
    })
}

#[derive(Clone)]
struct Frame {
    node: u64,
    next_mask: u32,                 // unexplored outgoing edges (bitmask)
    total_outdeg: usize,            // total outdegree of this node
    solved: usize,                  // # internal bifurcations crossed (exclude start)
    ref_species: BitVec<u32, Lsb0>, // cumulative species reference set for the path
    emitted: bool,                  // whether this path ending was already recorded
}

/// Iterative DFS that:
/// - Starts at start (a bifurcation), walks until an end_kmers node or depth limit.
/// - Initializes the species reference set from the start node and progressively refines it
///   using intersections involving the previous and current node at bifurcations.
/// - On reaching an end, records the accumulated species set with the path.
fn collect_paths_iterative_species_from_last_bifurcation(
    g: &Graph,
    start: u64,
    end_kmers: &HashSet<u64>,
    max_depth: usize,
    mask_k1: u64,
    min_path_edges: usize,
    out: &mut VariantGroups,
) {
    let start_mask = g.adj.get(start).unwrap_or(0);
    if start_mask == 0 {
        return;
    }
    let start_outdeg = start_mask.count_ones() as usize;

    let mut path: Vec<u64> = Vec::with_capacity(128);
    let mut on_path: HashSet<u64, RandomState> =
        HashSet::with_capacity_and_hasher(128, RandomState::new());

    let mut species_cache: HashMap<u64, BitVec<u32, Lsb0>, RandomState> =
        HashMap::with_capacity_and_hasher(128, RandomState::new());

    let mut get_species = |node: u64| {
        if let Some(bits) = species_cache.get(&node) {
            return bits.clone();
        }

        let bits = g
            .samples
            .get(node)
            .map(|bs| bs.to_bitvec())
            .unwrap_or_else(|| BitVec::repeat(false, g.n_species));

        species_cache.insert(node, bits.clone());
        bits
    };

    // Start species = species at start node
    let start_species = get_species(start);

    let start_frame = Frame {
        node: start,
        next_mask: start_mask,
        total_outdeg: start_outdeg,
        solved: 0,
        ref_species: start_species,
        emitted: false,
    };

    let mut stack: Vec<Frame> = Vec::with_capacity(256);
    stack.push(start_frame);
    path.push(start);
    on_path.insert(start);

    while let Some(mut fr) = stack.pop() {
        // Reached an end and path long enough: Emit path, materializing the species only now.
        if !fr.emitted && end_kmers.contains(&fr.node) && path.len() > min_path_edges {
            let species = fr.ref_species.clone();

            // In the node-based graph, the path is already explicit
            let full_nodes = path.clone();

            out.entry((start, fr.node)).or_default().push(PathRec {
                nodes: full_nodes,
                species,
            });

            fr.emitted = true;
        }

        if fr.next_mask != 0 {
            // Peel lowest set bit → deterministic successor order
            let b = fr.next_mask.trailing_zeros();
            fr.next_mask &= !(1 << b);

            // Put back the current frame with updated mask
            let cur = fr.node;
            stack.push(fr);

            // Compute successor directly from current node
            let next = ((cur << RECODE_BITS_PER_SYMBOL) | (b as u64)) & mask_k1;

            // Guard cycles cheaply
            if on_path.contains(&next) {
                continue;
            }

            // Depth control: increment when *leaving* an internal bifurcation (not the start)
            let parent = stack.last().unwrap();
            let cur_outdeg = parent.total_outdeg;
            let add = if parent.node != start && cur_outdeg > 1 {
                1
            } else {
                0
            };
            let next_solved = parent.solved + add;

            if next_solved > max_depth {
                continue;
            }

            // ---------- NEW SPECIES LOGIC ----------
            let mut next_ref_species = parent.ref_species.clone();

            if parent.node == start {
                // First step: species = S(start) ∩ S(next)
                let next_species = get_species(next);
                if next_ref_species.len() < next_species.len() {
                    next_ref_species.resize(next_species.len(), false);
                }
                next_ref_species &= next_species.as_bitslice();
            } else if parent.total_outdeg > 1 {
                // At a real bifurcation:
                //   species_new = current_ref ∩ S(parent) ∩ S(next)
                let parent_species = get_species(parent.node);
                if next_ref_species.len() < parent_species.len() {
                    next_ref_species.resize(parent_species.len(), false);
                }
                next_ref_species &= parent_species.as_bitslice();

                let next_species = get_species(next);
                if next_ref_species.len() < next_species.len() {
                    next_ref_species.resize(next_species.len(), false);
                }
                next_ref_species &= next_species.as_bitslice();
            } else {
                // Non-branching internal node:
                // keep current ref_species as-is for now.
            }
            // ---------- END NEW SPECIES LOGIC ----------

            // Prepare child frame
            let child_mask = g.adj.get(next).unwrap_or(0);
            let child_outdeg = child_mask.count_ones() as usize;
            let child = Frame {
                node: next,
                next_mask: child_mask,
                total_outdeg: child_outdeg,
                solved: next_solved,
                ref_species: next_ref_species,
                emitted: false,
            };

            stack.push(child);
            path.push(next);
            on_path.insert(next);
            continue;
        }

        // Backtrack
        if let Some(popped) = path.pop() {
            on_path.remove(&popped);
        }
    }
}

/// Per-group filtering:
///  (i) pick a single path length by max union species (break ties → drop group),
/// (ii) require ≥2 paths at that length,
/// (iii) require union species ≥ need.
fn filter_one_group_3steps(
    _g: &Graph,
    key: (u64, u64),
    paths: Vec<PathRec>,
    need: usize,
    n_species: usize,
) -> Option<((u64, u64), Vec<PathRec>)> {
    // (i) Choose path length by union species across paths of that length
    let mut unions_by_len: HashMap<usize, BitVec<u32, Lsb0>> = HashMap::new();
    for p in &paths {
        if p.nodes.len() < 2 {
            continue;
        }
        let len_edges = p.nodes.len() - 1;
        let u = unions_by_len
            .entry(len_edges)
            .or_insert_with(|| BitVec::repeat(false, n_species));
        if u.len() < p.species.len() {
            u.resize(p.species.len(), false);
        }
        *u |= p.species.as_bitslice();
    }
    if unions_by_len.is_empty() {
        return None;
    }

    // pick by (coverage desc, length asc). Drop ties on coverage.
    let mut items: Vec<(usize, usize)> = unions_by_len
        .iter()
        .map(|(len, bv)| (*len, bv.count_ones()))
        .collect();
    items.sort_unstable_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

    if items.len() >= 2 && items[0].1 == items[1].1 {
        return None; // ambiguous best length
    }
    let keep_len = items[0].0;

    // Keep only paths of the chosen length
    let kept_paths: Vec<PathRec> = paths
        .into_iter()
        .filter(|p| p.nodes.len() >= 2 && (p.nodes.len() - 1) == keep_len)
        .collect();

    // (ii) require ≥2 paths
    if kept_paths.len() < 2 {
        return None;
    }

    // (iii) union species across kept paths must be ≥ need
    let mut union_all: BitVec<u32, Lsb0> = BitVec::repeat(false, n_species);
    for p in &kept_paths {
        if union_all.len() < p.species.len() {
            union_all.resize(p.species.len(), false);
        }
        union_all |= p.species.as_bitslice();
    }
    if union_all.count_ones() < need {
        return None;
    }

    Some((key, kept_paths))
}

/// Score, sort, and keep non-overlapping (start,end) groups.
fn post_select_variant_groups(mut groups: VariantGroups) -> VariantGroups {
    let mut scored: Vec<ScoredGroup> = Vec::new();

    // Compute score (number of paths / path length)
    for (key, paths) in groups.drain() {
        let path_len = paths[0].nodes.len().saturating_sub(1); // all paths have the same length
        let n_paths = paths.len();
        let score = n_paths as f64 / path_len as f64;

        scored.push(((key, paths), score, path_len));
    }

    // Sort by score desc, then (start, end) for determinism
    scored.sort_unstable_by(|a, b| {
        // a.1 = score_a, b.1 = score_b
        b.1.partial_cmp(&a.1)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.0 .0 .0.cmp(&b.0 .0 .0)) // then by start node
            .then_with(|| a.0 .0 .1.cmp(&b.0 .0 .1)) // then by end node
    });

    let mut used_starts: HashSet<u64, RandomState> = HashSet::with_hasher(RandomState::new());
    let mut used_ends: HashSet<u64, RandomState> = HashSet::with_hasher(RandomState::new());
    let mut kept: VariantGroups = VariantGroups::default();

    for ((key, paths), _score, _len) in scored {
        let (s, e) = key;
        if used_starts.contains(&s) || used_ends.contains(&e) {
            continue;
        }
        used_starts.insert(s);
        used_ends.insert(e);
        kept.insert(key, paths);
    }

    kept
}
