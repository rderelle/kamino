#!/usr/bin/env bash
set -euo pipefail

INPUT_DIR="tests/data/input"
EXPECTED_DIR="tests/data/expected"
WORK_DIR="target/test_3diff"

mkdir -p "$WORK_DIR"

# Run from a clean working directory so outputs land there
(
  cd "$WORK_DIR"
  ../../target/release/kamino -i ../../"$INPUT_DIR"
)

# Compare the three expected output files with the program output
diff -u "${EXPECTED_DIR}/kamino_alignment.fas" "${WORK_DIR}/kamino_alignment.fas"
diff -u "${EXPECTED_DIR}/kamino_missing.tsv"    "${WORK_DIR}/kamino_missing.tsv"
diff -u "${EXPECTED_DIR}/kamino_partitions.tsv" "${WORK_DIR}/kamino_partitions.tsv"

echo "Golden test passed."
