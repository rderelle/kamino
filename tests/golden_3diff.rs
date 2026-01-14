use std::fs;
use std::path::Path;

use kamino_cli::{run_with_args, Args};

#[test]
fn golden_3diff() {
    let input_dir = Path::new("tests/data/input");
    let expected_dir = Path::new("tests/data/expected");
    let work_dir = Path::new("target/test_3diff_rs");

    // Make sure the directories we expect exist
    assert!(
        input_dir.exists(),
        "Input dir {:?} does not exist",
        input_dir
    );
    assert!(
        expected_dir.exists(),
        "Expected dir {:?} does not exist",
        expected_dir
    );

    // Create work directory
    fs::create_dir_all(work_dir).expect("Failed to create work dir");

    // Use an absolute path for input, but write outputs in work_dir
    let input_abs = input_dir
        .canonicalize()
        .expect("Failed to canonicalize input dir");
    let output_prefix = work_dir.join("kamino");

    let args = Args {
        input: Some(input_abs),
        input_file: None,
        k: None,
        min_freq: 0.85,
        depth: 4,
        output: output_prefix,
        constant: None,
        length_middle: None,
        mask: 5,
        threads: None,
        recode: kamino_cli::RecodeScheme::Dayhoff6,
        version: (),
    };

    run_with_args(args).expect("Failed to run kamino pipeline");

    // Helper to compare two files like `diff`
    fn compare_files(expected: &Path, got: &Path) {
        let expected_content = fs::read_to_string(expected).unwrap_or_else(|e| {
            panic!("Failed to read expected file {:?}: {}", expected, e);
        });
        let got_content = fs::read_to_string(got).unwrap_or_else(|e| {
            panic!("Failed to read output file {:?}: {}", got, e);
        });

        assert_eq!(
            expected_content, got_content,
            "File contents differ: expected {:?}, got {:?}",
            expected, got
        );
    }

    // Compare the 3 expected output files with the program output
    compare_files(
        &expected_dir.join("kamino_alignment.fas"),
        &work_dir.join("kamino_alignment.fas"),
    );
    compare_files(
        &expected_dir.join("kamino_missing.tsv"),
        &work_dir.join("kamino_missing.tsv"),
    );
    compare_files(
        &expected_dir.join("kamino_partitions.tsv"),
        &work_dir.join("kamino_partitions.tsv"),
    );

    println!("Golden test passed (Rust integration test).");
}
