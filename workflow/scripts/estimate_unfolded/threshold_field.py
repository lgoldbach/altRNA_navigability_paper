#!/usr/bin/env python3
"""
threshold_field.py

Read a 'genotype value' file (as written by diffusion_unfolded_field.py)
and pick a threshold so that a target fraction of genotypes are selected
as "unfolded" (values above the threshold).

Usage as a library
-------------------
    from threshold_field import load_field, select_unfolded

    genotypes, values = load_field("field.txt")
    unfolded_genotypes, threshold, achieved_fraction = select_unfolded(
        genotypes, values, target_fraction=0.87)

Command-line usage
-------------------
    python threshold_field.py field.txt --target-fraction 0.87 \\
        --output unfolded_genotypes.txt
"""

import sys
import time
import argparse
import numpy as np


def load_field(path):
    """Load a 'genotype value' file into (genotype_array, value_array).
    Tries pandas first (much faster for ~10^7-line files); falls back to
    a manual parser if pandas isn't installed."""
    t0 = time.time()
    try:
        import pandas as pd
        df = pd.read_csv(path, sep=" ", header=None, names=["genotype", "value"],
                          dtype={"genotype": str, "value": np.float64},
                          engine="c")
        genotypes = df["genotype"].to_numpy()
        values = df["value"].to_numpy(dtype=np.float64)
    except ImportError:
        genotypes_list = []
        values_list = []
        with open(path) as f:
            for line in f:
                parts = line.split()
                if not parts:
                    continue
                genotypes_list.append(parts[0])
                values_list.append(float(parts[1]))
        genotypes = np.array(genotypes_list)
        values = np.array(values_list, dtype=np.float64)

    print(f"Loaded {len(genotypes):,} genotypes from {path} "
          f"({time.time()-t0:.1f}s)", file=sys.stderr)
    return genotypes, values


def pick_threshold(values, target_fraction):
    """Return the threshold such that approximately target_fraction of
    values exceed it (values > threshold => 'unfolded')."""
    if not 0.0 <= target_fraction <= 1.0:
        raise ValueError("target_fraction must be in [0, 1]")
    quantile_level = 1.0 - target_fraction
    return float(np.quantile(values, quantile_level))


def select_unfolded(genotypes, values, target_fraction):
    """Return (unfolded_genotypes, threshold, achieved_fraction)."""
    threshold = pick_threshold(values, target_fraction)
    mask = values > threshold
    achieved_fraction = float(mask.mean())
    return genotypes[mask], threshold, achieved_fraction


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("field_file")
    ap.add_argument("--target-fraction", type=float, required=True,
                     help="Desired fraction of genotypes to mark unfolded, e.g. 0.87")
    ap.add_argument("--output", default=None,
                     help="Optional path to write the selected unfolded "
                          "genotype strings to, one per line.")
    args = ap.parse_args()

    genotypes, values = load_field(args.field_file)
    unfolded, threshold, achieved = select_unfolded(
        genotypes, values, args.target_fraction)

    print(f"threshold = {threshold:.6g}", file=sys.stderr)
    print(f"target fraction = {args.target_fraction:.4f}, "
          f"achieved fraction = {achieved:.4f} "
          f"({len(unfolded):,} genotypes)", file=sys.stderr)

    if args.output:
        with open(args.output, "w") as f:
            for g in unfolded:
                f.write(g + "\n")
        print(f"Wrote {len(unfolded):,} unfolded genotypes to {args.output}",
              file=sys.stderr)


if __name__ == "__main__":
    main()
