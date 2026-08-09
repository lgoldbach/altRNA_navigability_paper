#!/usr/bin/env python3
"""
genotypes_to_indices.py

Convert a file of genotype strings into a file of their integer indices
(0-indexed line numbers in a genotypes.txt file) -- the inverse of
indices_to_genotypes.py.

This does not assume genotypes.txt is sorted in any particular order
(e.g. it doesn't assume the seq_to_id base-4 encoding used elsewhere) --
it builds an explicit string -> line-number lookup from the file's actual
contents, so it works regardless of how genotypes.txt happens to be
ordered.

Implementation note: a plain Python dict mapping each of the ~1.6x10^7
genotype strings to its line number turns out to use too much memory for
a ~4GB machine (millions of individually-boxed Python string objects plus
hash-table overhead). Instead, this sorts a fixed-width numpy string
array once (genotypes.astype fixed to the observed sequence length, so
numpy can compare/sort them natively in C rather than falling back to
per-element Python object comparisons) and looks up each query via
vectorized binary search (np.searchsorted) -- both more memory-frugal and
faster than the dict approach at this scale.

Usage
------
    python genotypes_to_indices.py genotype_strings.txt genotypes.txt \\
        --output indices.txt
"""

import sys
import time
import argparse
import numpy as np
import pandas as pd


def load_genotypes_array(genotypes_path):
    """Load genotypes.txt into a fixed-width numpy string array indexed
    by (0-based) line number -- fixed width (rather than pandas' default
    Python-object string dtype) is what lets numpy sort/compare these
    natively in C."""
    df = pd.read_csv(genotypes_path, header=None, names=["genotype"],
                      dtype={"genotype": str}, engine="c")
    genotypes = df["genotype"].to_numpy(dtype=str)
    max_len = max((len(g) for g in genotypes), default=0)
    return genotypes.astype(f"<U{max_len}")


def build_sorted_lookup(genotypes_array):
    """Return (sorted_genotypes, sort_order) such that
    sorted_genotypes = genotypes_array[sort_order], enabling vectorized
    binary search for arbitrary query genotypes via np.searchsorted."""
    sort_order = np.argsort(genotypes_array, kind="stable")
    sorted_genotypes = genotypes_array[sort_order]
    return sorted_genotypes, sort_order


def lookup_indices(queries_array, sorted_genotypes, sort_order):
    """Vectorized lookup of each query's original (pre-sort) line index.
    Returns (indices, found_mask); indices is -1 wherever found_mask is
    False (the query genotype wasn't present in genotypes.txt at all).

    Note: if genotypes.txt contains duplicate lines, searchsorted returns
    the position of the first matching entry in sorted order, so a query
    matching a duplicated genotype resolves to whichever of its original
    line numbers happens to sort first -- genotypes.txt is not expected
    to have duplicates, but this is worth knowing if it does."""
    positions = np.searchsorted(sorted_genotypes, queries_array)
    positions_clipped = np.clip(positions, 0, len(sorted_genotypes) - 1)
    found_mask = sorted_genotypes[positions_clipped] == queries_array
    indices = np.full(len(queries_array), -1, dtype=np.int64)
    indices[found_mask] = sort_order[positions_clipped[found_mask]]
    return indices, found_mask


def convert(query_path, genotypes_path, output_path, chunk_size=500_000):
    t0 = time.time()

    genotypes = load_genotypes_array(genotypes_path)
    print(f"Loaded {len(genotypes):,} genotypes from {genotypes_path} "
          f"({time.time()-t0:.1f}s)", file=sys.stderr)

    sorted_genotypes, sort_order = build_sorted_lookup(genotypes)
    print(f"Built sorted lookup ({time.time()-t0:.1f}s)", file=sys.stderr)

    queries_df = pd.read_csv(query_path, header=None, names=["genotype"],
                              dtype={"genotype": str}, engine="c")
    max_len = max((len(g) for g in genotypes), default=0)
    queries = queries_df["genotype"].to_numpy(dtype=str).astype(f"<U{max_len}")
    print(f"Loaded {len(queries):,} query genotypes from {query_path} "
          f"({time.time()-t0:.1f}s)", file=sys.stderr)

    indices, found_mask = lookup_indices(queries, sorted_genotypes, sort_order)
    n_missing = int((~found_mask).sum())
    if n_missing:
        examples = queries[~found_mask][:5].tolist()
        print(f"WARNING: {n_missing:,} query genotypes were not found in "
              f"{genotypes_path}; written as -1; first few missing: "
              f"{examples}", file=sys.stderr)

    print(f"Looked up {len(indices):,} indices ({time.time()-t0:.1f}s)",
          file=sys.stderr)

    with open(output_path, "w") as f:
        for start in range(0, len(indices), chunk_size):
            chunk = indices[start:start + chunk_size]
            f.write("\n".join(map(str, chunk.tolist())) + "\n")

    print(f"Wrote {len(indices):,} indices to {output_path} "
          f"({time.time()-t0:.1f}s total)", file=sys.stderr)

    return n_missing


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("genotype_strings_file",
                     help="File of genotype strings to look up, one per line.")
    ap.add_argument("genotypes_file",
                     help="genotypes.txt, one genotype string per line; "
                          "line number (0-indexed) is the id.")
    ap.add_argument("--output", default="indices.txt")
    ap.add_argument("--strict", action="store_true",
                     help="Exit with an error if any query genotype is "
                          "not found in genotypes_file (default: warn and "
                          "write -1 for those entries).")
    args = ap.parse_args()

    n_missing = convert(args.genotype_strings_file, args.genotypes_file, args.output)
    if args.strict and n_missing:
        sys.exit(f"error: {n_missing} genotypes not found in "
                  f"{args.genotypes_file} (--strict mode)")


if __name__ == "__main__":
    main()
