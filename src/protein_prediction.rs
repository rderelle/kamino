use anyhow::{bail, Context, Result};
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use seq_io::fasta::{Reader, Record};
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;

use crate::io::SpeciesInput;

const DEFAULT_MIN_AA: usize = 90;
const MAX_ALLOWED_OVERLAP_NT: usize = 30;

const GC_ESTIMATE_MAX_VALID_BASES: usize = 300_000;

// Additional local suppression among starts inside the same stop-free segment.
// If two starts are very close, keep only the better one.
const MIN_START_SEPARATION_NT: usize = 12;

// keep at most this many valid starts per stop-free segment.
const START_FILTER_LIMIT: usize = 8;
const ALT_START_TOTAL_GAP: i32 = 6;
const ALT_START_MAX: usize = 4;
const SD_SCAN_UPSTREAM_MAX: usize = 20;
const SD_SCAN_UPSTREAM_MIN: usize = 4;
const ADAPTIVE_TRAIN_MIN_AA: usize = 120;
const ADAPTIVE_MIN_TRAIN_ORFS: usize = 60;
const ADAPTIVE_MAX_TRAIN_ORFS: usize = usize::MAX;
const ADAPTIVE_MAX_TRAIN_NT: usize = usize::MAX;
const ADAPTIVE_BACKGROUND_NT_PER_STRAND: usize = usize::MAX;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum StartCodon {
    Atg,
    Gtg,
    Ttg,
}

impl StartCodon {
    #[inline]
    fn bonus(self) -> i32 {
        match self {
            StartCodon::Atg => 8,
            StartCodon::Gtg => 4,
            StartCodon::Ttg => 1,
        }
    }

    #[inline]
    fn rank(self) -> i32 {
        match self {
            StartCodon::Atg => 0,
            StartCodon::Gtg => 1,
            StartCodon::Ttg => 1,
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct GcModel {
    // expected GC proportions at codon positions 1,2,3, scaled by 1000
    exp1_pm: u16,
    exp2_pm: u16,
    exp3_pm: u16,
}

impl GcModel {
    #[inline]
    fn from_genome_gc(genome_gc: f32) -> Self {
        // Empirical regressions from 20 annotated genomes
        let gc1 = (0.684392 * genome_gc + 0.233505).clamp(0.20, 0.90);
        let gc2 = (0.432433 * genome_gc + 0.183878).clamp(0.15, 0.75);
        let gc3 = (1.887194 * genome_gc - 0.394644).clamp(0.05, 0.98);

        Self {
            exp1_pm: (gc1 * 1000.0).round() as u16,
            exp2_pm: (gc2 * 1000.0).round() as u16,
            exp3_pm: (gc3 * 1000.0).round() as u16,
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct Orf {
    start: usize,
    end: usize,
    strand: i8,
    scan_start: usize,
    scan_end: usize,
    partial_5p: bool,
    partial_3p: bool,
    has_start_codon: bool,
    total_score: i32,
}

impl Orf {
    #[inline]
    fn len_nt(self) -> usize {
        self.end - self.start
    }

    #[inline]
    fn partial_class(self) -> u8 {
        match (self.partial_5p, self.partial_3p) {
            (false, false) => 0,
            (true, false) | (false, true) => 1,
            (true, true) => 2,
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct CandidateStart {
    pos: usize,
    start_codon: StartCodon,
    start_score: i32,
}

#[derive(Clone, Copy, Debug)]
struct ValidStart {
    cand: CandidateStart,
    len_nt: usize,
}

#[derive(Clone, Copy, Debug)]
struct RankedStart {
    vs: ValidStart,
    total: i32,
}

#[derive(Clone, Copy)]
struct SegmentEdges {
    touches_3p: bool,
}

#[derive(Clone, Debug)]
struct AdaptiveCodingModel {
    codon_score_pm: [i16; 64],
}

// Precomputed lookup tables for nucleotide encoding, complement,
// codon classification, and translation.
const INVALID_NT2: u8 = 255;
const NT2_LUT: [u8; 256] = build_nt2_lut();

const COMPLEMENT_LUT: [u8; 256] = build_complement_lut();

const CODON_TO_AA_LUT: [char; 64] = [
    'K', 'N', 'K', 'N', 'T', 'T', 'T', 'T', 'R', 'S', 'R', 'S', 'I', 'I', 'M', 'I', 'Q', 'H', 'Q',
    'H', 'P', 'P', 'P', 'P', 'R', 'R', 'R', 'R', 'L', 'L', 'L', 'L', 'E', 'D', 'E', 'D', 'A', 'A',
    'A', 'A', 'G', 'G', 'G', 'G', 'V', 'V', 'V', 'V', '*', 'Y', '*', 'Y', 'S', 'S', 'S', 'S', '*',
    'C', 'W', 'C', 'L', 'F', 'L', 'F',
];

const START_CODON_LUT: [Option<StartCodon>; 64] = build_start_codon_lut();
const IS_STOP_LUT: [bool; 64] = build_stop_lut();

const fn build_nt2_lut() -> [u8; 256] {
    let mut lut = [INVALID_NT2; 256];
    lut[b'A' as usize] = 0;
    lut[b'C' as usize] = 1;
    lut[b'G' as usize] = 2;
    lut[b'T' as usize] = 3;
    lut
}

const fn build_complement_lut() -> [u8; 256] {
    let mut lut = [b'N'; 256];
    lut[b'A' as usize] = b'T';
    lut[b'a' as usize] = b't';
    lut[b'T' as usize] = b'A';
    lut[b't' as usize] = b'a';
    lut[b'C' as usize] = b'G';
    lut[b'c' as usize] = b'g';
    lut[b'G' as usize] = b'C';
    lut[b'g' as usize] = b'c';
    lut
}

const fn build_start_codon_lut() -> [Option<StartCodon>; 64] {
    let mut lut = [None; 64];
    lut[14] = Some(StartCodon::Atg); // ATG
    lut[46] = Some(StartCodon::Gtg); // GTG
    lut[62] = Some(StartCodon::Ttg); // TTG
    lut
}

const fn build_stop_lut() -> [bool; 64] {
    let mut lut = [false; 64];
    lut[48] = true; // TAA
    lut[50] = true; // TAG
    lut[56] = true; // TGA
    lut
}

pub fn predict_proteomes(
    genomes: &[SpeciesInput],
    temp_dir: &Path,
    num_threads: usize,
) -> Result<Vec<SpeciesInput>> {
    let pool = ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .context("Failed to build Rayon thread pool for genome prediction")?;

    let predicted = pool.install(|| {
        genomes
            .par_iter()
            .map(|genome| -> Result<SpeciesInput> {
                let out_path = temp_dir.join(format!("{}.faa", genome.name));
                predict_one_genome(&genome.path, &out_path)?;
                Ok(SpeciesInput {
                    name: genome.name.clone(),
                    path: out_path,
                })
            })
            .collect::<Vec<Result<SpeciesInput>>>()
    });

    predicted.into_iter().collect()
}

// Main pipeline: parse arguments, scan each sequence, predict ORFs, and write translated proteins.
fn predict_one_genome(input_path: &Path, output_path: &Path) -> Result<()> {
    let min_orf_nt = DEFAULT_MIN_AA * 3;

    let (genome_gc, atgc_fraction) = estimate_input_genome_gc(input_path)?;
    if atgc_fraction < 0.5 {
        bail!(
            "the file {} does not appear to contain DNA sequences.",
            input_path.display()
        );
    }
    let gc_model = GcModel::from_genome_gc(genome_gc);

    let mut reader = open_fasta_reader(input_path)?;

    let out_file = File::create(output_path)
        .with_context(|| format!("Failed to create output file: {output_path:?}"))?;
    let mut out_writer = BufWriter::new(out_file);

    let mut protein_count = 0usize;
    let mut rc: Vec<u8> = Vec::new();
    let mut orfs: Vec<Orf> = Vec::new();

    while let Some(record_result) = reader.next() {
        let record = record_result.context("Failed to read FASTA record")?;

        let mut seq = record.seq().to_vec();
        seq.retain(|b| !b.is_ascii_whitespace());
        seq.make_ascii_uppercase();

        if seq.len() < 3 {
            continue;
        }

        reverse_complement_into(&seq, &mut rc);

        orfs.clear();
        find_orfs_on_strand(&seq, 1, min_orf_nt, &gc_model, &mut orfs);
        find_orfs_on_strand(&rc, -1, min_orf_nt, &gc_model, &mut orfs);

        if let Some(model) = train_adaptive_coding_model(&orfs, &seq, &rc) {
            apply_adaptive_coding_scores(&mut orfs, &seq, &rc, &model);
        }

        let accepted = greedy_non_overlapping(&mut orfs);

        for orf in accepted {
            let nt_slice = if orf.strand == 1 {
                &seq[orf.scan_start..orf.scan_end]
            } else {
                &rc[orf.scan_start..orf.scan_end]
            };

            let protein = translate_table11(nt_slice, orf.has_start_codon);

            protein_count += 1;
            write_one_protein_fasta(&mut out_writer, protein_count, &protein)?;
        }
    }

    out_writer
        .flush()
        .with_context(|| format!("Failed to flush output file: {output_path:?}"))?;

    Ok(())
}

fn open_fasta_reader(path: &Path) -> Result<Reader<Box<dyn Read>>> {
    let file = File::open(path).with_context(|| format!("Failed to open input file: {path:?}"))?;
    let buffered = BufReader::new(file);
    let reader: Box<dyn Read> = if path
        .extension()
        .and_then(|ext| ext.to_str())
        .is_some_and(|ext| ext.eq_ignore_ascii_case("gz"))
    {
        Box::new(flate2::read::MultiGzDecoder::new(buffered))
    } else {
        Box::new(buffered)
    };

    Ok(Reader::new(reader))
}

#[inline]
fn write_one_protein_fasta<W: Write>(writer: &mut W, idx: usize, protein: &str) -> Result<()> {
    writeln!(writer, ">protein_{idx}")
        .with_context(|| format!("Failed to write FASTA header for record {idx}"))?;

    for chunk in protein.as_bytes().chunks(60) {
        writer
            .write_all(chunk)
            .with_context(|| format!("Failed to write FASTA sequence for record {idx}"))?;
        writer
            .write_all(b"\n")
            .with_context(|| format!("Failed to write FASTA newline for record {idx}"))?;
    }

    Ok(())
}

// Estimate genome-wide GC content from the first valid bases of the input FASTA.
fn estimate_input_genome_gc(path: &Path) -> Result<(f32, f32)> {
    let mut reader = open_fasta_reader(path)?;

    let mut gc = 0usize;
    let mut valid = 0usize;
    let mut total = 0usize;

    while let Some(record_result) = reader.next() {
        if valid >= GC_ESTIMATE_MAX_VALID_BASES {
            break;
        }

        let record = record_result.context("Failed to read FASTA record during GC estimate")?;

        let mut seq = record.seq().to_vec();
        seq.retain(|b| !b.is_ascii_whitespace());
        
        for b in seq {
            if total >= GC_ESTIMATE_MAX_VALID_BASES {
                break;
            }

            total += 1;

            match b {
                b'G' | b'C' | b'g' | b'c' => {
                    gc += 1;
                    valid += 1;
                }
                b'A' | b'T' | b'a' | b't' => {
                    valid += 1;
                }
                _ => {}
            }
        }
    }

    let genome_gc = if valid == 0 {
        0.5
    } else {
        gc as f32 / valid as f32
    };

    let atgc_fraction = if total == 0 {
        0.0
    } else {
        valid as f32 / total as f32
    };

    Ok((genome_gc, atgc_fraction))
}

#[inline]
fn overlap_len(a_start: usize, a_end: usize, b_start: usize, b_end: usize) -> usize {
    let s = a_start.max(b_start);
    let e = a_end.min(b_end);
    e.saturating_sub(s)
}

// Keep the best-scoring ORFs while removing predictions that overlap too much.
#[inline]
fn greedy_non_overlapping(orfs: &mut [Orf]) -> Vec<Orf> {
    orfs.sort_unstable_by(|a, b| {
        a.partial_class()
            .cmp(&b.partial_class())
            .then_with(|| b.total_score.cmp(&a.total_score))
            .then_with(|| b.len_nt().cmp(&a.len_nt()))
            .then_with(|| a.start.cmp(&b.start))
    });

    let mut accepted: Vec<Orf> = Vec::with_capacity(orfs.len());
    for &cand in orfs.iter() {
        let mut ok = true;
        for a in &accepted {
            let ov = overlap_len(cand.start, cand.end, a.start, a.end);
            if ov > MAX_ALLOWED_OVERLAP_NT {
                ok = false;
                break;
            }
        }
        if ok {
            accepted.push(cand);
        }
    }
    accepted
}

// Scan one strand in all three frames, split it into stop-free segments, and collect ORF candidates.
fn find_orfs_on_strand(
    seq: &[u8],
    strand: i8,
    min_orf_nt: usize,
    gc_model: &GcModel,
    out: &mut Vec<Orf>,
) {
    let len = seq.len();
    let mut suppressed_starts: Vec<CandidateStart> = Vec::new();
    let mut valid_starts: Vec<ValidStart> = Vec::new();

    for frame in 0..3 {
        if frame >= len {
            continue;
        }

        let mut candidate_starts: Vec<CandidateStart> = Vec::new();
        let mut i = frame;
        while i + 3 <= len {
            let a = seq[i];
            let b = seq[i + 1];
            let c = seq[i + 2];
            let codon_has_n = has_n_triplet(a, b, c);

            if codon_has_n {
                emit_segment_orfs(
                    seq,
                    strand,
                    len,
                    frame,
                    i,
                    &candidate_starts,
                    min_orf_nt,
                    gc_model,
                    &mut suppressed_starts,
                    &mut valid_starts,
                    out,
                );
                candidate_starts.clear();
                i += 3;
                continue;
            }

            if let Some(sc) = start_codon_triplet(a, b, c) {
                let start_score = start_score(seq, i, sc);
                candidate_starts.push(CandidateStart {
                    pos: i,
                    start_codon: sc,
                    start_score,
                });
            }

            if is_stop_triplet(a, b, c) {
                emit_segment_orfs(
                    seq,
                    strand,
                    len,
                    frame,
                    i,
                    &candidate_starts,
                    min_orf_nt,
                    gc_model,
                    &mut suppressed_starts,
                    &mut valid_starts,
                    out,
                );
                candidate_starts.clear();
            }

            i += 3;
        }

        let frame_end = len - ((len - frame) % 3);
        emit_segment_orfs(
            seq,
            strand,
            len,
            frame,
            frame_end,
            &candidate_starts,
            min_orf_nt,
            gc_model,
            &mut suppressed_starts,
            &mut valid_starts,
            out,
        );
    }
}

// Evaluate one stop-free segment and emit the best ORF predictions from it.
#[inline]
#[allow(clippy::too_many_arguments)]
fn emit_segment_orfs(
    seq: &[u8],
    strand: i8,
    seq_len: usize,
    frame: usize,
    segment_end: usize,
    candidate_starts: &[CandidateStart],
    min_orf_nt: usize,
    gc_model: &GcModel,
    suppressed_starts: &mut Vec<CandidateStart>,
    valid_starts: &mut Vec<ValidStart>,
    out: &mut Vec<Orf>,
) {
    let frame_end = seq_len - ((seq_len - frame) % 3);
    let edges = SegmentEdges {
        touches_3p: segment_end == frame_end,
    };

    suppress_nearby_starts(candidate_starts, suppressed_starts);
    collect_valid_starts(suppressed_starts, segment_end, min_orf_nt, valid_starts);

    emit_complete_orfs(
        seq,
        strand,
        seq_len,
        segment_end,
        edges,
        gc_model,
        valid_starts,
        out,
    );
}

#[inline]
fn collect_valid_starts(
    starts: &[CandidateStart],
    segment_end: usize,
    min_orf_nt: usize,
    out: &mut Vec<ValidStart>,
) {
    out.clear();
    out.reserve(starts.len());

    for &s in starts {
        if s.pos >= segment_end {
            continue;
        }

        let len_nt = segment_end - s.pos;
        if len_nt >= min_orf_nt {
            out.push(ValidStart { cand: s, len_nt });
        }
    }
}

// Rank valid starts in a segment and keep only the best complete ORFs.
#[inline]
#[allow(clippy::too_many_arguments)]
fn emit_complete_orfs(
    seq: &[u8],
    strand: i8,
    seq_len: usize,
    segment_end: usize,
    edges: SegmentEdges,
    gc_model: &GcModel,
    valid_starts: &mut [ValidStart],
    out: &mut Vec<Orf>,
) -> bool {
    if valid_starts.is_empty() {
        return false;
    }

    let mut ranked: Vec<RankedStart> = Vec::with_capacity(valid_starts.len());
    for &vs in valid_starts.iter() {
        let nt_slice = &seq[vs.cand.pos..segment_end];
        let gc_score = gcpos_score_pm(nt_slice, gc_model);
        let total = total_orf_score(
            vs.len_nt,
            true,
            vs.cand.start_score,
            gc_score,
            false,
            edges.touches_3p,
        );
        ranked.push(RankedStart { vs, total });
    }

    ranked.sort_unstable_by(|a, b| {
        b.total
            .cmp(&a.total)
            .then_with(|| b.vs.cand.start_score.cmp(&a.vs.cand.start_score))
            .then_with(|| {
                a.vs.cand
                    .start_codon
                    .rank()
                    .cmp(&b.vs.cand.start_codon.rank())
            })
            .then_with(|| b.vs.len_nt.cmp(&a.vs.len_nt))
            .then_with(|| a.vs.cand.pos.cmp(&b.vs.cand.pos))
    });

    let top_total = ranked[0].total;

    for (idx, rs) in ranked.iter().enumerate() {
        if idx >= START_FILTER_LIMIT || idx >= ALT_START_MAX {
            break;
        }
        if idx > 0 && top_total - rs.total > ALT_START_TOTAL_GAP {
            break;
        }

        push_orf_with_total(
            strand,
            seq_len,
            rs.vs.cand.pos,
            segment_end,
            false,
            edges.touches_3p,
            true,
            rs.total,
            out,
        );
    }

    true
}

#[inline]
fn suppress_nearby_starts(starts: &[CandidateStart], kept: &mut Vec<CandidateStart>) {
    kept.clear();
    kept.reserve(starts.len());

    'outer: for &s in starts {
        for k in &mut *kept {
            let dist = s.pos.abs_diff(k.pos);
            if dist < MIN_START_SEPARATION_NT {
                let replace = s.start_score > k.start_score
                    || (s.start_score == k.start_score
                        && s.start_codon.rank() < k.start_codon.rank());

                if replace {
                    *k = s;
                }
                continue 'outer;
            }
        }
        kept.push(s);
    }
}

#[inline]
#[allow(clippy::too_many_arguments)]
fn push_orf_with_total(
    strand: i8,
    seq_len: usize,
    scan_start: usize,
    scan_end: usize,
    partial_5p: bool,
    partial_3p: bool,
    has_start_codon: bool,
    total_score: i32,
    out: &mut Vec<Orf>,
) {
    let (start, end) = if strand == 1 {
        (scan_start, scan_end)
    } else {
        (seq_len - scan_end, seq_len - scan_start)
    };

    out.push(Orf {
        start,
        end,
        strand,
        scan_start,
        scan_end,
        partial_5p,
        partial_3p,
        has_start_codon,
        total_score,
    });
}

#[inline]
fn start_score(seq: &[u8], start: usize, sc: StartCodon) -> i32 {
    sc.bonus() + sd_score(seq, start)
}

#[inline]
fn sd_score(seq: &[u8], start: usize) -> i32 {
    if start < SD_SCAN_UPSTREAM_MIN + 4 {
        return 0;
    }

    let lo = start.saturating_sub(SD_SCAN_UPSTREAM_MAX);
    let hi = start.saturating_sub(SD_SCAN_UPSTREAM_MIN);

    let mut best = 0;
    let mut i = lo;
    while i + 4 <= hi {
        let mut local = 0;
        if i + 5 <= hi {
            let mut m5 = 0;
            m5 += (seq[i] == b'A') as i32;
            m5 += (seq[i + 1] == b'G') as i32;
            m5 += (seq[i + 2] == b'G') as i32;
            m5 += (seq[i + 3] == b'A') as i32;
            m5 += (seq[i + 4] == b'G') as i32;

            local = match m5 {
                5 => 5,
                4 => 3,
                _ => 0,
            };
        }

        if local == 0 {
            let mut purines = 0;
            purines += is_purine(seq[i]) as i32;
            purines += is_purine(seq[i + 1]) as i32;
            purines += is_purine(seq[i + 2]) as i32;
            purines += is_purine(seq[i + 3]) as i32;
            local = match purines {
                4 => 1,
                3 => 0,
                _ => 0,
            };
        }

        if local > 0 {
            let motif_len = if local >= 3 { 5 } else { 4 };
            let spacer = start - (i + motif_len);
            local += match spacer {
                6..=10 => 2,
                5 | 11 => 1,
                _ => 0,
            };
        }

        best = best.max(local);
        i += 1;
    }

    best
}

#[inline]
fn total_orf_score(
    nt_len: usize,
    has_start_codon: bool,
    start_score: i32,
    gc_score: i32,
    partial_5p: bool,
    partial_3p: bool,
) -> i32 {
    let len_score = length_score(nt_len);
    let start_term = calibrated_start_term(start_score, nt_len);

    let mut score = 0i32;
    score += len_score;
    score += gc_score;
    score += start_term;

    if !has_start_codon {
        score -= 10;
    }
    if partial_5p {
        score -= 8;
    }
    if partial_3p {
        score -= 8;
    }

    score
}

#[inline]
fn calibrated_start_term(start_score: i32, nt_len: usize) -> i32 {
    if nt_len >= 600 {
        (start_score / 2).clamp(-6, 7)
    } else if nt_len >= 300 {
        ((2 * start_score) / 3).clamp(-8, 9)
    } else {
        start_score.clamp(-10, 10)
    }
}

#[inline]
fn length_score(nt_len: usize) -> i32 {
    let aa_len = nt_len / 3;

    if aa_len < 20 {
        return -20;
    }

    let x = aa_len as i32;
    (60 * x) / (x + 80)
}

#[inline]
fn gcpos_score_pm(seq: &[u8], gc_model: &GcModel) -> i32 {
    let codons = seq.len() / 3;
    if codons < 20 {
        return 0;
    }

    let usable = codons * 3;
    let s = &seq[..usable];

    let mut gc1 = 0usize;
    let mut gc2 = 0usize;
    let mut gc3 = 0usize;

    let mut i = 0usize;
    while i < usable {
        if is_gc(s[i]) {
            gc1 += 1;
        }
        if is_gc(s[i + 1]) {
            gc2 += 1;
        }
        if is_gc(s[i + 2]) {
            gc3 += 1;
        }
        i += 3;
    }

    let obs1_pm = ((gc1 * 1000) / codons) as i32;
    let obs2_pm = ((gc2 * 1000) / codons) as i32;
    let obs3_pm = ((gc3 * 1000) / codons) as i32;

    let d1 = (obs1_pm - gc_model.exp1_pm as i32).abs();
    let d2 = (obs2_pm - gc_model.exp2_pm as i32).abs();
    let d3 = (obs3_pm - gc_model.exp3_pm as i32).abs();

    let weighted_penalty = d1 + d2 + 2 * d3;

    let mut score = 90 - (weighted_penalty / 25);
    score += (codons.min(400) / 40) as i32;

    score
}

// Build a simple codon-usage model from high-confidence ORFs in the current sequence.
#[inline]
fn train_adaptive_coding_model(
    orfs: &[Orf],
    fwd: &[u8],
    rev: &[u8],
) -> Option<AdaptiveCodingModel> {
    let mut seed: Vec<Orf> = orfs
        .iter()
        .copied()
        .filter(|o| {
            !o.partial_5p
                && !o.partial_3p
                && o.has_start_codon
                && o.len_nt() / 3 >= ADAPTIVE_TRAIN_MIN_AA
        })
        .collect();
    if seed.len() < ADAPTIVE_MIN_TRAIN_ORFS {
        return None;
    }

    seed.sort_unstable_by(|a, b| b.total_score.cmp(&a.total_score));
    let keep = (seed.len() / 2)
        .clamp(ADAPTIVE_MIN_TRAIN_ORFS, ADAPTIVE_MAX_TRAIN_ORFS)
        .min(seed.len());
    seed.truncate(keep);

    let mut coding = [1usize; 64];
    let mut background = [1usize; 64];

    let mut trained_nt = 0usize;
    for orf in seed {
        let remaining_nt = ADAPTIVE_MAX_TRAIN_NT.saturating_sub(trained_nt);
        if remaining_nt == 0 {
            break;
        }
        let seq = if orf.strand == 1 { fwd } else { rev };
        let nt = &seq[orf.scan_start..orf.scan_end];
        let take = nt.len().min(remaining_nt);
        let usable = (take / 3) * 3;
        if usable >= 3 {
            accumulate_codon_counts(&nt[..usable], &mut coding);
            trained_nt += usable;
        }
    }

    let fwd_usable = (fwd.len().min(ADAPTIVE_BACKGROUND_NT_PER_STRAND) / 3) * 3;
    let rev_usable = (rev.len().min(ADAPTIVE_BACKGROUND_NT_PER_STRAND) / 3) * 3;
    if fwd_usable >= 3 {
        accumulate_codon_counts(&fwd[..fwd_usable], &mut background);
    }
    if rev_usable >= 3 {
        accumulate_codon_counts(&rev[..rev_usable], &mut background);
    }

    let coding_total: usize = coding.iter().sum();
    let background_total: usize = background.iter().sum();

    let mut codon_score_pm = [0i16; 64];
    for i in 0..64 {
        let p_coding = coding[i] as f64 / coding_total as f64;
        let p_bg = background[i] as f64 / background_total as f64;
        let lod2 = (p_coding / p_bg).log2();
        codon_score_pm[i] = (lod2 * 1000.0).clamp(-2200.0, 2200.0) as i16;
    }

    Some(AdaptiveCodingModel { codon_score_pm })
}

// Add the adaptive codon-usage score to each predicted ORF.
#[inline]
fn apply_adaptive_coding_scores(
    orfs: &mut [Orf],
    fwd: &[u8],
    rev: &[u8],
    model: &AdaptiveCodingModel,
) {
    for orf in orfs {
        let seq = if orf.strand == 1 { fwd } else { rev };
        let nt = &seq[orf.scan_start..orf.scan_end];
        let add = adaptive_coding_score(nt, model);
        orf.total_score += add;
    }
}

#[inline]
fn adaptive_coding_score(nt: &[u8], model: &AdaptiveCodingModel) -> i32 {
    let codons = nt.len() / 3;
    if codons < 35 {
        return 0;
    }

    let usable = (nt.len() / 3) * 3;
    let mut sum_pm = 0i32;
    let mut i = 0usize;
    while i < usable {
        if let Some(idx) = codon_index_fast(nt[i], nt[i + 1], nt[i + 2]) {
            sum_pm += model.codon_score_pm[idx] as i32;
        }
        i += 3;
    }

    let avg_pm = sum_pm / codons as i32;
    (avg_pm / 25).clamp(-18, 18)
}

#[inline]
fn accumulate_codon_counts(seq: &[u8], counts: &mut [usize; 64]) {
    let usable = (seq.len() / 3) * 3;
    let mut i = 0usize;
    while i < usable {
        if let Some(idx) = codon_index_fast(seq[i], seq[i + 1], seq[i + 2]) {
            counts[idx] += 1;
        }
        i += 3;
    }
}

#[inline]
fn codon_index_fast(a: u8, b: u8, c: u8) -> Option<usize> {
    let x = NT2_LUT[a as usize];
    let y = NT2_LUT[b as usize];
    let z = NT2_LUT[c as usize];
    if x == INVALID_NT2 || y == INVALID_NT2 || z == INVALID_NT2 {
        return None;
    }
    Some(((x << 4) | (y << 2) | z) as usize)
}

// Build the reverse-complemented version of the input sequence for scanning the opposite strand.
fn reverse_complement_into(seq: &[u8], rc: &mut Vec<u8>) {
    rc.clear();
    rc.resize(seq.len(), b'N');

    let len = seq.len();
    for (i, &b) in seq.iter().enumerate() {
        rc[len - 1 - i] = COMPLEMENT_LUT[b as usize];
    }
}

#[inline]
fn has_n_triplet(a: u8, b: u8, c: u8) -> bool {
    a == b'N' || b == b'N' || c == b'N'
}

#[inline]
fn is_stop_triplet(a: u8, b: u8, c: u8) -> bool {
    if let Some(idx) = codon_index_fast(a, b, c) {
        IS_STOP_LUT[idx]
    } else {
        false
    }
}

#[inline]
fn start_codon_triplet(a: u8, b: u8, c: u8) -> Option<StartCodon> {
    if let Some(idx) = codon_index_fast(a, b, c) {
        START_CODON_LUT[idx]
    } else {
        None
    }
}

#[inline]
fn is_purine(b: u8) -> bool {
    b == b'A' || b == b'G'
}

#[inline]
fn is_gc(b: u8) -> bool {
    b == b'G' || b == b'C'
}

// Translate a nucleotide ORF into a protein sequence using bacterial translation table 11.
fn translate_table11(nt: &[u8], force_initial_m: bool) -> String {
    let aa_len = nt.len() / 3;
    let mut protein = String::with_capacity(aa_len);

    let mut i = 0usize;
    if force_initial_m && nt.len() >= 3 {
        protein.push('M');
        i = 3;
    }

    while i + 3 <= nt.len() {
        let aa = codon_to_aa_table11(nt[i], nt[i + 1], nt[i + 2]);
        protein.push(aa);
        i += 3;
    }

    protein
}

#[inline]
fn codon_to_aa_table11(a: u8, b: u8, c: u8) -> char {
    if let Some(idx) = codon_index_fast(a, b, c) {
        CODON_TO_AA_LUT[idx]
    } else {
        'X'
    }
}
