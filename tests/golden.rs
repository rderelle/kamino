//! Golden integration tests for checked-in fixture datasets.
//!
//! Each case runs the public library API on fixture FASTA files and compares the
//! deterministic output files against the checked-in expected results.

use std::fs;
use std::path::Path;

use kamino_cli::{run_with_args, Args, RecodeScheme};

mod common;

struct GoldenCase {
    name: &'static str,
    input_dir: Option<&'static str>,
    input_file: Option<&'static str>,
    expected_dir: &'static str,
    work_dir: &'static str,
    genomes: bool,
    length_middle: usize,
    recode: RecodeScheme,
    nj: bool,
}

fn run_golden_case(case: GoldenCase) {
    let expected_dir = Path::new(case.expected_dir);
    let work_dir = Path::new(case.work_dir);

    if let Some(input_dir) = case.input_dir {
        common::assert_fixture_dirs(Path::new(input_dir), expected_dir);
    }
    if let Some(input_file) = case.input_file {
        assert!(
            Path::new(input_file).exists(),
            "Input file {:?} does not exist",
            input_file
        );
    }

    if work_dir.exists() {
        fs::remove_dir_all(work_dir).unwrap_or_else(|e| {
            panic!(
                "Failed to remove previous work dir {:?} for golden case {}: {}",
                work_dir, case.name, e
            );
        });
    }
    fs::create_dir_all(work_dir).unwrap_or_else(|e| {
        panic!(
            "Failed to create work dir {:?} for golden case {}: {}",
            work_dir, case.name, e
        );
    });

    // Use absolute paths for inputs, but write outputs in work_dir.
    let input_abs = case.input_dir.map(|input_dir| {
        Path::new(input_dir).canonicalize().unwrap_or_else(|e| {
            panic!(
                "Failed to canonicalize input dir {:?} for golden case {}: {}",
                input_dir, case.name, e
            );
        })
    });
    let input_file_abs = case.input_file.map(|input_file| {
        Path::new(input_file).canonicalize().unwrap_or_else(|e| {
            panic!(
                "Failed to canonicalize input file {:?} for golden case {}: {}",
                input_file, case.name, e
            );
        })
    });
    let output_prefix = work_dir.join("kamino");

    let args = Args {
        input: input_abs,
        input_file: input_file_abs,
        genomes: case.genomes,
        k: None,
        min_freq: 0.85,
        output: output_prefix,
        constant: None,
        length_middle: case.length_middle,
        mask: 5,
        threads: 1,
        recode: case.recode,
        nj: case.nj,
    };

    run_with_args(args)
        .unwrap_or_else(|e| panic!("Failed to run kamino pipeline for {}: {}", case.name, e));
    common::compare_test_outputs(expected_dir, work_dir);

    println!("Golden test passed: {} (Rust integration test).", case.name);
}

#[test]
fn golden_3diff() {
    run_golden_case(GoldenCase {
        name: "3diff",
        input_dir: Some("tests/data/test_3diff/input"),
        input_file: None,
        expected_dir: "tests/data/test_3diff/expected",
        work_dir: "target/test_3diff_rs",
        genomes: false,
        length_middle: 35,
        recode: RecodeScheme::Dayhoff6,
        nj: true,
    });
}

#[test]
fn golden_3diff_input_file_kgb6() {
    run_golden_case(GoldenCase {
        name: "3diff_input_file_kgb6",
        input_dir: None,
        input_file: Some("tests/data/test_3diff/input.tsv"),
        expected_dir: "tests/data/test_3diff/expected",
        work_dir: "target/test_3diff_input_file_kgb6_rs",
        genomes: false,
        length_middle: 35,
        recode: RecodeScheme::KGB6,
        nj: true,
    });
}

#[test]
fn golden_genomes() {
    run_golden_case(GoldenCase {
        name: "genomes",
        input_dir: Some("tests/data/test_genomes/input"),
        input_file: None,
        expected_dir: "tests/data/test_genomes/expected",
        work_dir: "target/test_genomes_rs",
        genomes: true,
        length_middle: 50,
        recode: RecodeScheme::SR6,
        nj: false,
    });
}

#[test]
fn proteome_output_is_identical_across_thread_counts_and_repeated_runs() {
    let input = Path::new("tests/data/test_3diff/input")
        .canonicalize()
        .unwrap();
    let temp = tempfile::tempdir().unwrap();
    let output_files = [
        "kamino_alignment.fas",
        "kamino_missing.tsv",
        "kamino_partitions.tsv",
    ];
    let mut outputs = Vec::new();

    for (run, threads) in [1, 4, 4].into_iter().enumerate() {
        let output_dir = temp.path().join(format!("run_{run}"));
        fs::create_dir(&output_dir).unwrap();
        run_with_args(Args {
            input: Some(input.clone()),
            input_file: None,
            genomes: false,
            k: None,
            min_freq: 0.85,
            output: output_dir.join("kamino"),
            constant: None,
            length_middle: 35,
            mask: 5,
            threads,
            recode: RecodeScheme::Dayhoff6,
            nj: false,
        })
        .unwrap();
        outputs.push(
            output_files
                .iter()
                .map(|name| fs::read(output_dir.join(name)).unwrap())
                .collect::<Vec<_>>(),
        );
    }

    assert_eq!(
        outputs[0], outputs[1],
        "1-thread and 4-thread output differs"
    );
    assert_eq!(outputs[1], outputs[2], "repeated 4-thread output differs");
}
