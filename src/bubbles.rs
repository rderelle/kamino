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

    // Shared helpers
    #[inline]
    fn bump_epoch(epoch: &mut u32, marks: &mut [u32]) {
        *epoch = epoch.wrapping_add(1);
        if *epoch == 0 {
            marks.fill(0);
            *epoch = 1;
        }
    }

    #[inline]
    fn union_count_until(
        marks: &mut [u32],
        epoch: u32,
        colors: NodeColorSlice,
        need: usize,
        mut covered: usize,
    ) -> usize {
        for sid in colors.iter_ones() {
            let s = sid;
            if unsafe { *marks.get_unchecked(s) } != epoch {
                unsafe { *marks.get_unchecked_mut(s) = epoch };
                covered += 1;
                if covered >= need {
                    break;
                }
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
        // START nodes (outdeg ≥ 2) — match END logic:
        // union of species across outgoing paths
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

                    succs.clear();
                    scratch.clear();

                    let mut m = mask;
                    while m != 0 {
                        let b = m.trailing_zeros() as u64;
                        m &= m - 1;
                        let v = ((u << RECODE_BITS_PER_SYMBOL) | b) & main.k1_mask;
                        succs.push(v);
                    }

                    let mut sum = 0usize;
                    for &v in succs.iter() {
                        let cnt = main.samples.get(v).map(|bs| bs.count_ones()).unwrap_or(0);
                        sum += cnt;
                        scratch.push((v, cnt));
                    }

                    if sum < need {
                        return None;
                    }

                    scratch.sort_unstable_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

                    bump_epoch(epoch, marks.as_mut_slice());
                    let mut covered = 0usize;
                    for (v, _c) in scratch.iter().copied() {
                        if let Some(bs) = main.samples.get(v) {
                            covered =
                                union_count_until(marks.as_mut_slice(), *epoch, bs, need, covered);
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
        // END nodes (indeg ≥ 2) — unchanged logic:
        // union of species across incoming paths
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

                    let tail = (v & sym_mask) as u32;
                    let base = v >> RECODE_BITS_PER_SYMBOL;
                    let shift_bits_u128 = shift_bits as u128;
                    if shift_bits_u128 >= 128 {
                        return None;
                    }

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

                    let mut sum = 0usize;
                    for &u in preds.iter() {
                        let cnt = main.samples.get(u).map(|bs| bs.count_ones()).unwrap_or(0);
                        sum += cnt;
                        scratch.push((u, cnt));
                    }
                    if sum < need {
                        return None;
                    }

                    scratch.sort_unstable_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

                    bump_epoch(epoch, marks.as_mut_slice());
                    let mut covered = 0usize;
                    for (u, _c) in scratch.iter().copied() {
                        if let Some(bs) = main.samples.get(u) {
                            covered =
                                union_count_until(marks.as_mut_slice(), *epoch, bs, need, covered);
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
