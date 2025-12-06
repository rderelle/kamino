use anyhow::{Context, Result};
use dashmap::DashMap;
use flate2::read::MultiGzDecoder;
use hashbrown::HashMap;
use seq_io::fasta::{Reader as FastaReader, Record};
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};
use std::hash::BuildHasherDefault;
use rustc_hash::FxHasher;

use crate::graph::{AdjTable, Graph, NodeColorTable, SpeciesSet};
use crate::recode::{recode_byte, RecodeScheme, RECODE_BITS_PER_SYMBOL};

use rayon::prelude::*;
use rayon::ThreadPoolBuilder;

// Type alias for faster hashing
type FastHashMap<K, V> = HashMap<K, V, BuildHasherDefault<FxHasher>>;
type FastDashMap<K, V> = DashMap<K, V, BuildHasherDefault<FxHasher>>;

// ------------------------------
// File I/O helpers
// ------------------------------

pub fn collect_fasta_files(input_dir: &Path) -> anyhow::Result<Vec<PathBuf>> {
    fn is_supported_ext(p: &Path) -> bool {
        let ext = p.extension().and_then(|e| e.to_str()).unwrap_or("");
        if ext.eq_ignore_ascii_case("gz") {
            if let Some(stem) = p.file_stem() {
                let stem_path = Path::new(stem);
                let inner = stem_path.extension().and_then(|e| e.to_str()).unwrap_or("");
                matches_ignore_ascii(inner, &["fa", "fas", "fasta", "faa"])
            } else {
                false
            }
        } else {
            matches_ignore_ascii(ext, &["fa", "fas", "fasta", "faa"])
        }
    }

    #[inline]
    fn matches_ignore_ascii(ext: &str, set: &[&str]) -> bool {
        set.iter().any(|&e| ext.eq_ignore_ascii_case(e))
    }

    let mut files: Vec<PathBuf> = std::fs::read_dir(input_dir)?
        .filter_map(|e| e.ok())
        .filter(|e| e.file_type().map(|t| t.is_file()).unwrap_or(false))
        .map(|e| e.path())
        .filter(|p| is_supported_ext(p))
        .collect();

    files.sort_unstable_by(|a, b| a.file_name().unwrap().cmp(b.file_name().unwrap()));

    Ok(files)
}

pub(crate) fn species_name_from_path(p: &std::path::Path) -> String {
    let mut stem = p
        .file_name()
        .and_then(|x| x.to_str())
        .unwrap_or("unnamed")
        .to_string();
    if stem.ends_with(".gz") {
        stem.truncate(stem.len() - 3);
    }
    for ext in [".fa", ".fasta", ".fas", ".faa"] {
        if stem.to_ascii_lowercase().ends_with(ext) {
            let n = stem.len() - ext.len();
            stem.truncate(n);
            break;
        }
    }
    stem
}

pub(crate) fn open_fasta(path: &Path) -> Result<Box<dyn Read>> {
    let f = File::open(path).with_context(|| format!("open {:?}", path))?;
    // Larger buffer for better I/O throughput
    let buffered = BufReader::with_capacity(512 * 1024, f);
    let r: Box<dyn Read> = if let Some(s) = path.to_str() {
        if s.ends_with(".gz") {
            Box::new(MultiGzDecoder::new(buffered))
        } else {
            Box::new(buffered)
        }
    } else {
        Box::new(buffered)
    };
    Ok(r)
}

// ------------------------------
// Graph builder
// ------------------------------

pub fn build_graph_from_dir(
    dir: &Path,
    k: usize,
    main: &mut Graph,
    num_threads: usize,
) -> Result<()> {
    let files = collect_fasta_files(dir)?;
    eprintln!("proteome files: {}", files.len());

    main.init_species_len(files.len());

    // Pre-register species sequentially for deterministic IDs
    let mut species_ids: Vec<u16> = Vec::with_capacity(files.len());
    for path in &files {
        let sid = main.register_species(species_name_from_path(path)) as u16;
        species_ids.push(sid);
    }

    let sym_bits = RECODE_BITS_PER_SYMBOL;
    let k_mask = main.k_mask;
    let k1_mask = main.k1_mask;
    let recode_scheme = main.recode_scheme;
    let sym_mask: u64 = if sym_bits >= 64 {
        u64::MAX
    } else {
        (1u64 << sym_bits) - 1
    };

    // Single-pass kmer streaming
    #[inline]
    #[allow(clippy::too_many_arguments)]
    fn stream_kmers<R: Read>(
        mut reader: FastaReader<R>,
        k: usize,
        sym_bits: u32,
        k_mask: u64,
        k1_mask: u64,
        sym_mask: u64,
        mut on_edge: impl FnMut(u64, u32, u64),
        scheme: RecodeScheme,
    ) -> Result<()> {
        while let Some(rec) = reader.next() {
            let rec = rec?;
            let seq = rec.seq();

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
                    let prefix = (roll >> sym_bits) & k1_mask;
                    let sym = (roll & sym_mask) as u32;
                    let suffix = ((prefix << sym_bits) | (sym as u64)) & k1_mask;
                    on_edge(prefix, sym, suffix);
                }
            }
        }
        Ok(())
    }

    let pool = ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .expect("Failed to build local Rayon thread pool");

    let adj: FastDashMap<u64, u32> = FastDashMap::default();
    let samples: FastDashMap<u64, SpeciesSet> = FastDashMap::default();

    pool.install(|| -> Result<()> {
        files
            .par_iter()
            .zip(species_ids.par_iter().copied())
            .try_for_each(|(path, sid)| -> Result<()> {
                let rdr = open_fasta(path)?;
                let reader = FastaReader::new(rdr);

                // Thread-local adjacency tracking (with_capacity needs further optimisation)
                let mut local_adj: FastHashMap<u64, u32> = FastHashMap::with_capacity_and_hasher(
                    1_000_000,
                    BuildHasherDefault::<FxHasher>::default()
                );

                // Single pass: build adjacency masks
                stream_kmers(
                    reader,
                    k,
                    sym_bits,
                    k_mask,
                    k1_mask,
                    sym_mask,
                    |u, sym, _v| {
                        let bit = 1u32 << sym;
                        *local_adj.entry(u).or_insert(0) |= bit;
                    },
                    recode_scheme,
                )?;

                let mut kept = 0usize;
                let mut dropped = 0usize;

                // Pre-allocate based on expected filtering
                let estimated = local_adj.len();  // 95%+ pass rate
                let mut adj_updates: Vec<(u64, u32)> = Vec::with_capacity(estimated);
                let mut node_updates: Vec<u64> = Vec::with_capacity(estimated * 2);

                // Filter: keep only nodes with out-degree == 1
                for (&u, &mask) in &local_adj {
                    if mask == 0 {
                        continue;
                    }

                    let od = mask.count_ones();
                    if od != 1 {
                        dropped += od as usize;
                        continue;
                    }

                    let sym = mask.trailing_zeros();
                    let bit = 1u32 << sym;
                    let v = ((u << sym_bits) | (sym as u64)) & k1_mask;

                    adj_updates.push((u, bit));
                    node_updates.push(u);
                    node_updates.push(v);
                    kept += 1;
                }

                // Drop local_adj early to free memory
                drop(local_adj);

                eprintln!(
                    "[{}] kept={} dropped={}",
                    path.file_name().and_then(|x| x.to_str()).unwrap_or("?"),
                    kept,
                    dropped
                );

                // Batch merge adjacency
                for (node, mask) in adj_updates {
                    adj.entry(node)
                        .and_modify(|m| *m |= mask)
                        .or_insert(mask);
                }

                // Sort and dedup nodes
                node_updates.sort_unstable();
                node_updates.dedup();

                // Batch update samples
                for node in node_updates {
                    samples.entry(node)
                        .and_modify(|species_set| {
                            if let Err(pos) = species_set.binary_search(&sid) {
                                species_set.insert(pos, sid);
                            }
                        })
                        .or_insert_with(|| {
                            let mut set = SpeciesSet::new();
                            set.push(sid);
                            set
                        });
                }

                Ok(())
            })
    })?;

    main.adj = AdjTable::from_iter(adj);
    main.samples = NodeColorTable::from_iter(samples, main.n_species);

    Ok(())
}

pub fn print_graph_size(g: &Graph) {
    let nodes = g.adj.len();
    let edges: usize = g
        .adj
        .iter()
        .map(|(_, mask)| mask.count_ones() as usize)
        .sum();
    eprintln!("graph size: nodes={} edges={}", nodes, edges);
}
