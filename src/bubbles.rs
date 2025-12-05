use hashbrown::HashSet;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;

use crate::graph::{Graph, NodeColorSlice};
use crate::recode::{RECODE_ALPHABET_SIZE, RECODE_BITS_PER_SYMBOL};

pub fn find_bubble_endpoints(
    main: &Graph,
    ratio: f32,
    num_threads: usize,
) -> (HashSet<u64>, HashSet<u64>) {
    let n_species = main.n_species.max(1);
    let need = ((ratio * n_species as f32).ceil() as usize).min(n_species);

    let mut sources: Vec<u64> = main.adj.nodes().to_vec();
    sources.sort_unstable();

    let sym_mask = (1u64 << RECODE_BITS_PER_SYMBOL) - 1;
    let shift_symbols = main.k.saturating_sub(2);
    let shift_bits = RECODE_BITS_PER_SYMBOL.saturating_mul(shift_symbols as u32);

    // ----------------------------
    // Small helpers
    // ----------------------------

    #[inline]
    fn bump_epoch(epoch: &mut u32, marks: &mut [u32]) {
        *epoch = epoch.wrapping_add(1);
        if *epoch == 0 {
            marks.fill(0);
            *epoch = 1;
        }
    }

    /// |colors_u ∩ colors_v| where each is a sorted set of species IDs
    #[inline]
    fn intersection_count(colors_u: NodeColorSlice, colors_v: NodeColorSlice) -> usize {
        let mut it_u = colors_u.iter_ones();
        let mut it_v = colors_v.iter_ones();

        let mut a = it_u.next();
        let mut b = it_v.next();
        let mut count = 0usize;

        while let (Some(x), Some(y)) = (a, b) {
            if x == y {
                count += 1;
                a = it_u.next();
                b = it_v.next();
            } else if x < y {
                a = it_u.next();
            } else {
                b = it_v.next();
            }
        }

        count
    }

    /// Union of intersections:
    /// covered species in ⋃ (colors_u ∩ colors_v_i), written into `marks`.
    #[inline]
    fn union_intersection_until(
        marks: &mut [u32],
        epoch: u32,
        colors_u: NodeColorSlice,
        colors_v: NodeColorSlice,
        need: usize,
        mut covered: usize,
    ) -> usize {
        // two-pointer intersection on the fly; mark species in `marks`
        let mut it_u = colors_u.iter_ones();
        let mut it_v = colors_v.iter_ones();

        let mut a = it_u.next();
        let mut b = it_v.next();

        while let (Some(x), Some(y)) = (a, b) {
            if x == y {
                let s = x; // species ID (usize)
                unsafe {
                    if *marks.get_unchecked(s) != epoch {
                        *marks.get_unchecked_mut(s) = epoch;
                        covered += 1;
                        if covered >= need {
                            break;
                        }
                    }
                }
                a = it_u.next();
                b = it_v.next();
            } else if x < y {
                a = it_u.next();
            } else {
                b = it_v.next();
            }
        }

        covered
    }

    // Build a local Rayon thread pool with the requested number of threads
    let n_threads = num_threads.max(1);
    let pool = ThreadPoolBuilder::new()
        .num_threads(n_threads)
        .build()
        .expect("Failed to build Rayon thread pool in find_bubble_endpoints");

    // Run the parallel parts inside this pool
    pool.install(|| {
        // ============================================================
        // START nodes (outdeg ≥ 2)
        // Use union of intersections:
        //   S(u) = ⋃_v (colors(u) ∩ colors(v))
        // ============================================================

        let start_vec: Vec<u64> = sources
            .par_iter()
            .map_init(
                || {
                    (
                        vec![0u32; n_species],
                        Vec::with_capacity(RECODE_ALPHABET_SIZE),
                        Vec::with_capacity(RECODE_ALPHABET_SIZE),
                        1u32,
                    )
                },
                |state, &u| {
                    let (marks, succs, scratch, epoch) = state;

                    let mask = main.adj.get(u).unwrap_or(0);
                    if mask.count_ones() < 2 {
                        return None;
                    }

                    // colours of the start node itself
                    let colors_u = match main.samples.get(u) {
                        Some(bs) => bs,
                        None => return None,
                    };

                    succs.clear();
                    scratch.clear();

                    // collect successors of u
                    let mut m = mask;
                    while m != 0 {
                        let b = m.trailing_zeros() as u64;
                        m &= m - 1;
                        let v = ((u << RECODE_BITS_PER_SYMBOL) | b) & main.k1_mask;
                        succs.push(v);
                    }

                    // quick upper-bound check using intersection sizes
                    let mut sum = 0usize;
                    for &v in succs.iter() {
                        let cnt = if let Some(colors_v) = main.samples.get(v) {
                            intersection_count(colors_u, colors_v)
                        } else {
                            0
                        };
                        sum += cnt;
                        scratch.push((v, cnt));
                    }

                    if sum < need {
                        return None;
                    }

                    // sort successors by intersection size, then by node id
                    scratch.sort_unstable_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

                    bump_epoch(epoch, marks.as_mut_slice());
                    let mut covered = 0usize;

                    // S(u) = ⋃_v (colors(u) ∩ colors(v))
                    for (v, _c) in scratch.iter().copied() {
                        if let Some(colors_v) = main.samples.get(v) {
                            covered = union_intersection_until(
                                marks.as_mut_slice(),
                                *epoch,
                                colors_u,
                                colors_v,
                                need,
                                covered,
                            );
                            if covered >= need {
                                return Some(u);
                            }
                        }
                    }

                    None
                },
            )
            .filter_map(|x| x)
            .collect();

        let start_kmers: HashSet<u64> = start_vec.into_iter().collect();

        // =========================================
        // END nodes (indeg ≥ 2)
        // Symmetric logic using predecessors:
        //   S(v) = ⋃_u (colors(v) ∩ colors(u))
        // =========================================

        let mut vs: Vec<u64> = main.samples.nodes().to_vec();
        vs.sort_unstable();
        vs.dedup();

        let end_vec: Vec<u64> = vs
            .par_iter()
            .map_init(
                || {
                    (
                        vec![0u32; n_species],
                        Vec::with_capacity(RECODE_ALPHABET_SIZE),
                        Vec::with_capacity(RECODE_ALPHABET_SIZE),
                        1u32,
                    )
                },
                |state, &v| {
                    let (marks, preds, scratch, epoch) = state;
                    preds.clear();
                    scratch.clear();

                    let colors_v = match main.samples.get(v) {
                        Some(bs) => bs,
                        None => return None,
                    };

                    let tail = (v & sym_mask) as u32;
                    let base = v >> RECODE_BITS_PER_SYMBOL;
                    let shift_bits_u128 = shift_bits as u128;
                    if shift_bits_u128 >= 128 {
                        return None;
                    }

                    // enumerate predecessors of v
                    for head in 0..(RECODE_ALPHABET_SIZE as u32) {
                        let candidate = (((head as u128) << shift_bits_u128) | base as u128)
                            & (main.k1_mask as u128);
                        let u = candidate as u64;
                        if let Some(mask) = main.adj.get(u) {
                            if (mask & (1u32 << tail)) != 0 {
                                preds.push(u);
                            }
                        }
                    }

                    preds.sort_unstable();
                    preds.dedup();

                    if preds.len() < 2 {
                        return None;
                    }

                    // quick upper-bound check using intersection sizes
                    let mut sum = 0usize;
                    for &u in preds.iter() {
                        let cnt = if let Some(colors_u) = main.samples.get(u) {
                            intersection_count(colors_u, colors_v)
                        } else {
                            0
                        };
                        sum += cnt;
                        scratch.push((u, cnt));
                    }

                    if sum < need {
                        return None;
                    }

                    // sort predecessors by intersection size, then by node id
                    scratch.sort_unstable_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

                    bump_epoch(epoch, marks.as_mut_slice());
                    let mut covered = 0usize;

                    // S(v) = ⋃_u (colors(v) ∩ colors(u))
                    for (u, _c) in scratch.iter().copied() {
                        if let Some(colors_u) = main.samples.get(u) {
                            covered = union_intersection_until(
                                marks.as_mut_slice(),
                                *epoch,
                                colors_v,
                                colors_u,
                                need,
                                covered,
                            );
                            if covered >= need {
                                return Some(v);
                            }
                        }
                    }

                    None
                },
            )
            .filter_map(|x| x)
            .collect();

        let end_kmers: HashSet<u64> = end_vec.into_iter().collect();

        (start_kmers, end_kmers)
    })
}
