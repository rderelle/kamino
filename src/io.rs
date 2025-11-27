use anyhow::{Context, Result};
use dashmap::DashMap;
use flate2::read::MultiGzDecoder;
use hashbrown::HashMap;
use seq_io::fasta::{Reader as FastaReader, Record};
use std::fs::File;
use std::io::Read;
use std::path::{Path, PathBuf};

use crate::graph::{AdjTable, Graph, NodeColorTable, SpeciesSet};
use crate::recode::recode_byte;

use rayon::prelude::*;
use rayon::ThreadPoolBuilder;

// ------------------------------
// File I/O helpers
// ------------------------------

/// Return FASTA/proteome files from input_dir, sorted alphabetically by filename.
/// Supported extensions (case-insensitive): .fa, .fas, .fasta, .faa, .fna and their .gz variants.
/// Sorting is lexicographic on the final file name (case-sensitive to be fully deterministic).
pub fn collect_fasta_files(input_dir: &Path) -> anyhow::Result<Vec<PathBuf>> {
    fn is_supported_ext(p: &Path) -> bool {
        // Accept: .fa, .fasta, .fas, .faa (any case)
        // Also: .fa.gz, .fasta.gz, .fas.gz, .faa.gz
        let ext = p.extension().and_then(|e| e.to_str()).unwrap_or("");
        if ext.eq_ignore_ascii_case("gz") {
            // Inspect the "inner" extension from the file_stem
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

    // Deterministic order: lexicographic by file name (no ties expected)
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
    let r: Box<dyn Read> = if let Some(s) = path.to_str() {
        if s.ends_with(".gz") {
            Box::new(MultiGzDecoder::new(f))
        } else {
            Box::new(f)
        }
    } else {
        Box::new(f)
    };
    Ok(r)
}

// ------------------------------
// Build graph
// ------------------------------

pub fn build_graph_from_dir(
    dir: &Path,
    k: usize,
    main: &mut Graph,
    num_threads: usize, // user-defined number of Rayon threads
) -> Result<()> {
    let files = collect_fasta_files(dir)?;
    eprintln!("proteome files: {}", files.len());

    main.init_species_len(files.len());

    // Pre-register species sequentially to keep IDs deterministic.
    let mut species_ids: Vec<u16> = Vec::with_capacity(files.len());
    for path in &files {
        let sid = main.register_species(species_name_from_path(path)) as u16;
        species_ids.push(sid);
    }

    let sym_bits = main.sym_bits;
    let k_mask = main.k_mask;
    let k1_mask = main.k1_mask;
    let sym_mask: u64 = if sym_bits >= 64 {
        u64::MAX
    } else {
        (1u64 << sym_bits) - 1
    };

    // Small helper to iterate kmers from a FASTA reader body without allocation.
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
    ) -> Result<()> {
        while let Some(rec) = reader.next() {
            let rec = rec?;
            let seq = rec.seq();

            let mut roll: u64 = 0;
            let mut have: usize = 0;

            for &b in seq {
                let a = recode_byte(b);
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

    // Build a local Rayon pool with the requested number of threads.
    let pool = ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .expect("Failed to build local Rayon thread pool");

    let adj: DashMap<u64, u32> = DashMap::new();
    let samples: DashMap<u64, SpeciesSet> = DashMap::new();

    // Parallel over proteome files inside this pool and update the global structures on the fly.
    pool.install(|| -> Result<()> {
        files
            .par_iter()
            .zip(species_ids.par_iter().copied())
            .try_for_each(|(path, sid)| -> Result<()> {
                let rdr = open_fasta(path)?;
                let reader = FastaReader::new(rdr);

                #[derive(Default)]
                struct NodeDeg {
                    out_mask: u32, // OR of outgoing symbol bits
                    indeg: u8,     // saturated at 2: 0,1,2 (2 = ">=2")
                }
                
                // with_capacity is enough for bacteria but not for eukaryotes
                let mut node_deg: HashMap<u64, NodeDeg> = HashMap::with_capacity(1_000_000);
                
                // First pass: fill node_deg with per-node out_mask and indegree (saturated).
                stream_kmers(
                    reader,
                    k,
                    sym_bits,
                    k_mask,
                    k1_mask,
                    sym_mask,
                    |u, sym, v| {
                        let bit = 1u32 << sym;

                        // Update out_mask for prefix node u
                        let entry_u = node_deg.entry(u).or_default();
                        if (entry_u.out_mask & bit) == 0 {
                            entry_u.out_mask |= bit;

                            // Only when this (u, sym) is new do we increment indegree for v
                            let entry_v = node_deg.entry(v).or_default();
                            if entry_v.indeg < 2 {
                                entry_v.indeg += 1; // 0 -> 1 -> 2 (2 == ">=2")
                            }
                        }
                    },
                )?;

                let mut kept = 0usize;
                let mut dropped = 0usize;

                let update_node = |node: u64| {
                    let mut entry = samples.entry(node).or_default();
                    let species = entry.value_mut();

                    match species.binary_search(&sid) {
                        Ok(_) => {
                            // already present
                        }
                        Err(pos) => {
                            species.insert(pos, sid);
                        }
                    }
                };

                // Second pass: filter edges based on od == 1 && indeg(v) == 1
                for (u, info) in node_deg.iter() {
                    let mask = info.out_mask;
                    if mask == 0 {
                        continue;
                    }

                    let od = mask.count_ones();

                    let mut bits = mask;
                    while bits != 0 {
                        let sym = bits.trailing_zeros();
                        bits &= bits - 1; // clear lowest set bit

                        let v = (((*u) << sym_bits) | (sym as u64)) & k1_mask;
                        let id = node_deg.get(&v).map(|nd| nd.indeg).unwrap_or(0);

                        if od == 1 && id == 1 {
                            let bit = 1u32 << sym;
                            adj.entry(*u)
                                .and_modify(|m| *m |= bit)
                                .or_insert(bit);

                            update_node(*u);
                            update_node(v);

                            kept += 1;
                        } else {
                            dropped += 1;
                        }
                    }
                }

                eprintln!(
                    "[{}] kept={} dropped={}",
                    path.file_name().and_then(|x| x.to_str()).unwrap_or("?"),
                    kept,
                    dropped
                );

                Ok(())
            })
    })?; // propagate any I/O / parsing errors

    main.adj = AdjTable::from_iter(adj, main.alphabet_size);
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

#[allow(dead_code)]
pub fn assert_adj_masks_within_alphabet(g: &Graph) {
    // Ensures we never set adjacency bits >= alphabet size.
    let max_bit = g.alphabet_size as u32;
    for (_, mask) in g.adj.iter() {
        if mask >> max_bit != 0 {
            panic!(
                "Adjacency has bits set >= alphabet_size ({}).",
                g.alphabet_size
            );
        }
    }
}
