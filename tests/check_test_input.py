#!/usr/bin/env python3

import os
import subprocess
import sys
import tempfile

import pandas as pd


def main():
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    script_path = os.path.join(repo_root, "peak_finder_batch.py")
    input_path = os.path.join(repo_root, "tests", "fixtures", "test_input.tsv")
    expected_path = os.path.join(repo_root, "tests", "fixtures", "test_input.expected.tsv")

    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = os.path.join(tmpdir, "test_output.tsv")
        subprocess.run(
            [sys.executable, script_path, "-f", input_path, "-o", output_path],
            check=True,
        )

        observed = pd.read_csv(output_path, sep="\t", dtype=str)
        expected = pd.read_csv(expected_path, sep="\t", dtype=str)

    if not observed.equals(expected):
        merged = observed.merge(
            expected,
            on=["trait", "chromosome", "bp_location", "rs_id", "pvalue"],
            how="outer",
            suffixes=("_obs", "_exp"),
        )
        diffs = merged[
            merged["isTopAssociation_obs"].fillna("") != merged["isTopAssociation_exp"].fillna("")
        ]
        raise AssertionError(
            "test_input regression failed.\n"
            + diffs.to_string(index=False)
        )

    print("test_input regression check passed.")


if __name__ == "__main__":
    main()
