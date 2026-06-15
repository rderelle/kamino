//! Input discovery and FASTA opening helpers.
//!
//! Inputs can come from a directory of FASTA files or from a two-column TSV table
//! mapping isolate names to FASTA paths. Files ending in `.gz` are decompressed on
//! the fly by [`open_fasta`].
use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

// ------------------------------
// File I/O helpers
// ------------------------------

pub fn collect_fasta_files(input_dir: &Path) -> anyhow::Result<Vec<PathBuf>> {
    // Accept both plain FASTA-like extensions and gzip-compressed versions such as
    // `.faa.gz`; sorting by filename keeps directory input deterministic.
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
/// One named species/proteome input consumed by the analysis pipeline.
pub struct SpeciesInput {
    /// Human-readable species or isolate name used in output FASTA headers.
    pub name: String,
    /// FASTA path for this species; may be gzip-compressed.
    pub path: PathBuf,
}

/// Build species inputs from all supported FASTA files in a directory.
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

/// Build species inputs from a two-column `name<TAB>path` table.
pub fn collect_species_inputs_from_table(table_path: &Path) -> anyhow::Result<Vec<SpeciesInput>> {
    // Relative paths in the table are resolved relative to the table file, not the
    // process working directory, which makes manifests portable.
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
    // Strip compression and FASTA extensions while preserving the rest of the file
    // name as the default isolate label.
    let mut stem = p
        .file_name()
        .and_then(|x| x.to_str())
        .unwrap_or("unnamed")
        .to_string();
    if stem.ends_with(".gz") {
        stem.truncate(stem.len() - 3);
    }
    for ext in [".fa", ".fasta", ".fas", ".faa", ".fna"] {
        if stem.to_ascii_lowercase().ends_with(ext) {
            let n = stem.len() - ext.len();
            stem.truncate(n);
            break;
        }
    }
    stem
}

pub(crate) fn open_fasta(path: &Path) -> Result<Box<dyn Read>> {
    // Return a boxed reader so callers can treat compressed and uncompressed FASTA
    // files identically.
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
