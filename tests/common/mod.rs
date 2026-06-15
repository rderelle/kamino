use std::fs;
use std::path::Path;

const DETERMINISTIC_OUTPUT_FILES: &[&str] = &[
    "kamino_alignment.fas",
    "kamino_missing.tsv",
    "kamino_partitions.tsv",
];

pub fn assert_fixture_dirs(input_dir: &Path, expected_dir: &Path) {
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
}

pub fn compare_test_outputs(expected_dir: &Path, work_dir: &Path) {
    for file_name in DETERMINISTIC_OUTPUT_FILES {
        compare_files(&expected_dir.join(file_name), &work_dir.join(file_name));
    }
}

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
