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
    input_dir: &'static str,
    expected_dir: &'static str,
    work_dir: &'static str,
    genomes: bool,
    length_middle: usize,
    recode: RecodeScheme,
    nj: bool,
}

fn run_golden_case(case: GoldenCase) {
    let input_dir = Path::new(case.input_dir);
    let expected_dir = Path::new(case.expected_dir);
    let work_dir = Path::new(case.work_dir);

    common::assert_fixture_dirs(input_dir, expected_dir);

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

    // Use an absolute path for input, but write outputs in work_dir.
    let input_abs = input_dir.canonicalize().unwrap_or_else(|e| {
        panic!(
            "Failed to canonicalize input dir {:?} for golden case {}: {}",
            input_dir, case.name, e
        );
    });
    let output_prefix = work_dir.join("kamino");

    let args = Args {
        input: Some(input_abs),
        input_file: None,
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
        input_dir: "tests/data/test_3diff/input",
        expected_dir: "tests/data/test_3diff/expected",
        work_dir: "target/test_3diff_rs",
        genomes: false,
        length_middle: 35,
        recode: RecodeScheme::Dayhoff6,
        nj: true,
    });
}

#[test]
fn golden_genomes() {
    run_golden_case(GoldenCase {
        name: "genomes",
        input_dir: "tests/data/test_genomes/input",
        expected_dir: "tests/data/test_genomes/expected",
        work_dir: "target/test_genomes_rs",
        genomes: true,
        length_middle: 50,
        recode: RecodeScheme::SR6,
        nj: false,
    });
}
