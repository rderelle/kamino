#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/.."

BIN="target/release/kamino"
DETERMINISTIC_OUTPUT_FILES=(
  "kamino_alignment.fas"
  "kamino_missing.tsv"
  "kamino_partitions.tsv"
)

run_case() {
  local name="$1"
  local input_dir="$2"
  local expected_dir="$3"
  local work_dir="$4"
  shift 4

  rm -rf "$work_dir"
  mkdir -p "$work_dir"

  # Run from a clean working directory so outputs land there.
  (
    cd "$work_dir"
    "../../${BIN}" "$@" -i "../../${input_dir}"
  )

  # Compare the deterministic expected output files with the program output.
  for file_name in "${DETERMINISTIC_OUTPUT_FILES[@]}"; do
    diff -u "${expected_dir}/${file_name}" "${work_dir}/${file_name}"
  done

  echo "Golden test passed: ${name}."
}

run_table_case() {
  local name="$1"
  local input_file="$2"
  local expected_dir="$3"
  local work_dir="$4"
  shift 4

  rm -rf "$work_dir"
  mkdir -p "$work_dir"

  # Run from a clean working directory so outputs land there.
  (
    cd "$work_dir"
    "../../${BIN}" "$@" -I "../../${input_file}"
  )

  # Compare the deterministic expected output files with the program output.
  for file_name in "${DETERMINISTIC_OUTPUT_FILES[@]}"; do
    diff -u "${expected_dir}/${file_name}" "${work_dir}/${file_name}"
  done

  echo "Golden test passed: ${name}."
}

run_case \
  "3diff" \
  "tests/data/test_3diff/input" \
  "tests/data/test_3diff/expected" \
  "target/test_3diff" \
  --length-middle 35 \
  --recode dayhoff6 \
  --nj

run_table_case \
  "3diff-input-file-kgb6" \
  "tests/data/test_3diff/input.tsv" \
  "tests/data/test_3diff/expected" \
  "target/test_3diff_input_file_kgb6" \
  --length-middle 35 \
  --recode kgb6 \
  --nj

run_case \
  "genomes" \
  "tests/data/test_genomes/input" \
  "tests/data/test_genomes/expected" \
  "target/test_genomes" \
  --genomes \
  --length-middle 50 \
  --recode sr6
