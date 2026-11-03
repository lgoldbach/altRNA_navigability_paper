#!/usr/bin/env python3
"""
indices_to_genotypes.py

Convert a file of integer genotype indices (0-indexed line numbers into a
genotypes.txt file -- index 0 is genotypes.txt's first line, and so on)
into a file of the corresponding genotype strings, one per line.

Usage
------
    python indices_to_genotypes.py unfolded_indices.txt genotypes.txt \\
        --output unfolded_genotypes.txt

The output file is directly usable as the true_unfolded_file input to
calibrate_diffusion.py.
"""

import sys
import time
import argparse
import numpy as np
import pandas as pd


def load_genotypes_array(genotypes_path):
    """Load genotypes.txt into a numpy string array indexed by (0-based)
    line number, for fast vectorized lookup."""
    df = pd.read_csv(genotypes_path, header=None, names=["genotype"],
                      dtype={"genotype": str}, engine="c")
    return df["genotype"].to_numpy()


def load_indices(indices_path):
    """Load a file of one non-negative integer per line into a numpy
    int64 array."""
    df = pd.read_csv(indices_path, header=None, names=["index"],
                      dtype={"index": np.int64}, engine="c")
    return df["index"].to_numpy()


def convert(indices_path, genotypes_path, output_path, chunk_size=500_000):
    t0 = time.time()

    genotypes = load_genotypes_array(genotypes_path)
    print(f"Loaded {len(genotypes):,} genotypes from {genotypes_path} "
          f"({time.time()-t0:.1f}s)", file=sys.stderr)

    indices = load_indices(indices_path)
    print(f"Loaded {len(indices):,} indices from {indices_path} "
          f"({time.time()-t0:.1f}s)", file=sys.stderr)

    if len(indices) == 0:
        raise ValueError(f"no indices found in {indices_path}")
    if indices.min() < 0:
        raise ValueError(f"found a negative index ({indices.min()}); "
                          f"indices must start at 0")
    max_idx = int(indices.max())
    if max_idx >= len(genotypes):
        raise ValueError(
            f"index {max_idx} is out of range for {genotypes_path}, which "
            f"has {len(genotypes):,} lines (valid indices: 0..{len(genotypes)-1})")

    selected = genotypes[indices]  # vectorized lookup, no Python-level loop
    print(f"Looked up {len(selected):,} genotype strings "
          f"({time.time()-t0:.1f}s)", file=sys.stderr)

    with open(output_path, "w") as f:
        for start in range(0, len(selected), chunk_size):
            chunk = selected[start:start + chunk_size]
            f.write("\n".join(chunk.tolist()) + "\n")

    print(f"Wrote {len(selected):,} genotype strings to {output_path} "
          f"({time.time()-t0:.1f}s total)", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("indices_file",
                     help="File of integer genotype indices, one per line "
                          "(0-indexed).")
    ap.add_argument("genotypes_file",
                     help="genotypes.txt, one genotype string per line; "
                          "line number (0-indexed) is the id referenced "
                          "by indices_file.")
    ap.add_argument("--output", default="unfolded_genotypes.txt")
    args = ap.parse_args()
    convert(args.indices_file, args.genotypes_file, args.output)


if __name__ == "__main__":
    main()
