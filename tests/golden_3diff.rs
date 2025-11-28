use std::fs;
use std::path::Path;
use std::process::Command;

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

    // Path to the compiled kamino binary (provided by Cargo for tests)
    let kamino_exe = env!("CARGO_BIN_EXE_kamino");

    // Use an absolute path for input, but run from work_dir
    let input_abs = input_dir
        .canonicalize()
        .expect("Failed to canonicalize input dir");

    let status = Command::new(kamino_exe)
        .current_dir(work_dir)
        .arg("-i")
        .arg(input_abs.to_str().expect("Input path not valid UTF-8"))
        .status()
        .expect("Failed to run kamino");

    assert!(
        status.success(),
        "kamino exited with non-zero status: {:?}",
        status.code()
    );

    // Helper to compare two files like `diff`
    fn compare_files(expected: &Path, got: &Path) {
        let expected_content =
            fs::read_to_string(expected).unwrap_or_else(|e| {
                panic!("Failed to read expected file {:?}: {}", expected, e);
            });
        let got_content =
            fs::read_to_string(got).unwrap_or_else(|e| {
                panic!("Failed to read output file {:?}: {}", got, e);
            });

        assert_eq!(
            expected_content, got_content,
            "File contents differ: expected {:?}, got {:?}",
            expected, got
        );
    }

    // Compare the three expected output files with the program output
    compare_files(
        &expected_dir.join("kamino_alignment.fas"),
        &work_dir.join("kamino_alignment.fas"),
    );
    compare_files(
        &expec
