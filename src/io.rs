use anyhow::{Context, Result};
use dashmap::DashMap;
use flate2::read::MultiGzDecoder;
use rustc_hash::FxHasher;
use seq_io::fasta::{Reader as FastaReader, Record};
use std::fs::File;
use std::hash::BuildHasherDefault;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

use crate::graph::{AdjTable, Graph, NodeColorTable, SpeciesSet};
use crate::proba_filter::{
    protein_passes_counts, stream_recoded_edges, stream_recoded_kmers, AtomicCountMinSketch,
    BloomFilter, BLOOM_BITS, CMS_DEPTH, CMS_WIDTH,
};
use crate::recode::RECODE_BITS_PER_SYMBOL;

use rayon::prelude::*;
use rayon::ThreadPoolBuilder;

// Type alias for faster hashing
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
                matches_ignore_ascii(inner, &["fa", "fas", "fasta", "faa", "fna"])
            } else {
                false
            }
        } else {
            matches_ignore_ascii(ext, &["fa", "fas", "fasta", "faa", "fna"])
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

#[derive(Debug, Clone)]
pub struct SpeciesInput {
    pub name: String,
    pub path: PathBuf,
}

pub fn collect_species_inputs_from_dir(input_dir: &Path) -> anyhow::Result<Vec<SpeciesInput>> {
    let files = collect_fasta_files(input_dir)?;
    let mut seen = hashbrown::HashSet::new();
    let mut inputs: Vec<SpeciesInput> = Vec::with_capacity(files.len());

    for path in files {
        let name = species_name_from_path(&path);
        if !seen.insert(name.clone()) {
            anyhow::bail!(
                "isolate name {:?} appears more than once in input directory.",
                name
            );
        }
        inputs.push(SpeciesInput { name, path });
    }

    inputs.sort_unstable_by(|a, b| a.name.cmp(&b.name).then_with(|| a.path.cmp(&b.path)));
    Ok(inputs)
}

pub fn collect_species_inputs_from_table(table_path: &Path) -> anyhow::Result<Vec<SpeciesInput>> {
    let file = File::open(table_path).with_context(|| format!("open {:?}", table_path))?;
    let reader = BufReader::new(file);
    let base_dir = table_path.parent().unwrap_or_else(|| Path::new("."));

    let mut seen = hashbrown::HashSet::new();
    let mut inputs: Vec<SpeciesInput> = Vec::new();

    for (line_idx, line) in reader.lines().enumerate() {
        let line = line.with_context(|| format!("read {:?}", table_path))?;
        let trimmed = line.trim_end();
        if trimmed.is_empty() {
            continue;
        }

        let mut parts = trimmed.split('\t');
        let Some(raw_name) = parts.next() else {
            anyhow::bail!(
                "Line {} in {:?} does not contain a species name.",
                line_idx + 1,
                table_path
            );
        };
        let Some(raw_path) = parts.next() else {
            anyhow::bail!(
                "Line {} in {:?} is missing a proteome path.",
                line_idx + 1,
                table_path
            );
        };
        if parts.next().is_some() {
            anyhow::bail!(
                "Line {} in {:?} has more than two tab-delimited columns.",
                line_idx + 1,
                table_path
            );
        }

        let name = raw_name.trim();
        let path_str = raw_path.trim();
        if name.is_empty() || path_str.is_empty() {
            anyhow::bail!(
                "Line {} in {:?} must contain non-empty isolate name and path.",
                line_idx + 1,
                table_path
            );
        }

        if !seen.insert(name.to_string()) {
            anyhow::bail!(
                "isolate name {:?} appears more than once in {:?}.",
                name,
                table_path
            );
        }

        let mut path = PathBuf::from(path_str);
        if path.is_relative() {
            path = base_dir.join(path);
        }
        let metadata =
            std::fs::metadata(&path).with_context(|| format!("read metadata for {:?}", path))?;
        anyhow::ensure!(
            metadata.is_file(),
            "Proteome path {:?} (line {} in {:?}) is not a file.",
            path,
            line_idx + 1,
            table_path
        );

        inputs.push(SpeciesInput {
            name: name.to_string(),
            path,
        });
    }

    inputs.sort_unstable_by(|a, b| a.name.cmp(&b.name).then_with(|| a.path.cmp(&b.path)));
    Ok(inputs)
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

pub fn build_graph_from_inputs(
    inputs: &[SpeciesInput],
    k: usize,
    min_freq: f32,
    length_middle: usize,
    main: &mut Graph,
    num_threads: usize,
) -> Result<()> {
    eprintln!("proteome files: {}", inputs.len());

    // Register species with deterministic IDs
    main.init_species_len(inputs.len());
    let mut species_ids: Vec<u16> = Vec::with_capacity(inputs.len());
    for input in inputs {
        let sid = main.register_species(input.name.clone()) as u16;
        species_ids.push(sid);
    }

    let sym_bits = RECODE_BITS_PER_SYMBOL;
    let k1_mask = main.k1_mask;
    let recode_scheme = main.recode_scheme;

    let k1 = k.saturating_sub(1);
    let min_needed = (min_freq * inputs.len().max(1) as f32).ceil() as u32;
    let min_needed_cms = min_needed;

    let pool = ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .expect("Failed to build local Rayon thread pool");

    // Pass 1: Build Count-Min Sketch (species frequency of each (k-1)-mer)
    let cms_atomic = AtomicCountMinSketch::new(CMS_WIDTH, CMS_DEPTH);
    pool.install(|| -> Result<()> {
        inputs.par_iter().try_for_each(|input| -> Result<()> {
            let rdr = open_fasta(&input.path)?;
            let mut reader = FastaReader::new(rdr);
            let mut species_bloom = BloomFilter::new(BLOOM_BITS);

            while let Some(rec) = reader.next() {
                let rec = rec?;
                let seq = rec.seq();
                stream_recoded_kmers(seq, k1, sym_bits, k1_mask, recode_scheme, |kmer| {
                    if !species_bloom.test_and_set(kmer) {
                        cms_atomic.increment(kmer);
                    }
                });
            }
            Ok(())
        })
    })?;
    let cms = cms_atomic.snapshot();

    // Shared structures for final graph
    let adj: FastDashMap<u64, u32> = FastDashMap::default();
    let samples: FastDashMap<u64, SpeciesSet> = FastDashMap::default();

    // Pass 2: Main parallel processing (filter + edges + species tracking)
    pool.install(|| -> Result<()> {
        inputs
            .par_iter()
            .zip(species_ids.par_iter().copied())
            .try_for_each(|(input, sid)| -> Result<()> {
                let rdr = open_fasta(&input.path)?;
                let mut reader = FastaReader::new(rdr);

                // Collect all candidate edges from passing proteins
                let mut transitions: Vec<(u64, u32)> = Vec::with_capacity(1_000_000);

                let mut passed = 0usize;
                let mut filtered = 0usize;

                while let Some(rec) = reader.next() {
                    let rec = rec?;
                    let seq = rec.seq();

                    // Filter protein using CMS
                    if !protein_passes_counts(
                        seq,
                        k1,
                        sym_bits,
                        k1_mask,
                        recode_scheme,
                        min_needed_cms,
                        length_middle,
                        |kmer| cms.estimate(kmer),
                    ) {
                        filtered += 1;
                        continue;
                    }

                    passed += 1;
                    // Stream edges into temporary vector
                    stream_recoded_edges(
                        seq,
                        k1,
                        sym_bits,
                        k1_mask,
                        recode_scheme,
                        |u, sym, _v| {
                            let bit = 1u32 << sym;
                            transitions.push((u, bit));
                        },
                    );
                }

                // Optimisation 1: Aggregate edges with sort + linear scan
                let mut kept = 0usize;
                let mut dropped = 0usize;

                transitions.sort_unstable_by_key(|(u, _)| *u);

                let mut adj_updates: Vec<(u64, u32)> = Vec::with_capacity(1_000_000);
                let mut node_updates: Vec<u64> = Vec::with_capacity(2_000_000);

                let mut i = 0;
                while i < transitions.len() {
                    let u = transitions[i].0;
                    let mut mask = 0u32;
                    while i < transitions.len() && transitions[i].0 == u {
                        mask |= transitions[i].1;
                        i += 1;
                    }
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

                drop(transitions);

                let proteins_total = passed + filtered;
                let kmers_total = kept + dropped;
                eprintln!(
                    "[{}] proteins={}/{} kmers={}/{}",
                    input.name, passed, proteins_total, kept, kmers_total
                );

                // Merge edges into shared adjacency table
                for (node, mask) in adj_updates {
                    adj.entry(node).and_modify(|m| *m |= mask).or_insert(mask);
                }

                // Optimisation 2: Collect species per node (push + final sort)
                node_updates.sort_unstable();
                node_updates.dedup();

                for node in node_updates {
                    samples
                        .entry(node)
                        .and_modify(|species_set| species_set.push(sid))
                        .or_insert_with(|| {
                            let mut set = SpeciesSet::new();
                            set.push(sid);
                            set
                        });
                }

                Ok(())
            })
    })?;

    // Sort all SpeciesSets
    samples.iter_mut().for_each(|mut entry| {
        entry.value_mut().sort_unstable();
    });

    // Build final graph structures
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
