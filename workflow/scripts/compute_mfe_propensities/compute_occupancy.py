#!/usr/bin/env python3
"""
compute_occupancy.py

Compute the propensity-weighted occupancy of every phenotype for every
genotype, i.e. for genotype g and phenotype S in g's permissible set:

    occupancy(g, S) = prop(S) / sum(prop(S') for S' in permissible_set(g))

which approximates "fraction of time genotype g spends in structure S",
using the mfe propensities as a Boltzmann-like weight (see prop_Si in the
manuscript's Bradley-Terry model).

Inputs
------
permissible_sets.txt : one line per unique phenotype:
    <phenotype_dotbracket> <gid_1> <gid_2> ... <gid_k>
    (gid_i are integers referencing line numbers of genotypes.txt)

genotypes.txt : one genotype string per line. Only the line COUNT is used
    (to size the accumulator array); indexing convention (0- or 1-based)
    does not matter here since gids are only ever used as array indices
    and echoed back unchanged.

mfe_propensities.csv : one line per phenotype, space-separated:
    <phenotype_dotbracket> <propensity_float>

Output
------
phenotype_occupancy.txt : one line per phenotype that has a known
    propensity (a subset of permissible_sets.txt's lines, same order):
    <phenotype_dotbracket> <gid_1>:<occ_1> <gid_2>:<occ_2> ...

    Genotypes not listed on a phenotype's line are implicitly occupancy 0
    for that phenotype (the phenotype is not in their permissible set, or
    the genotype's permissible set contains no scored phenotypes at all).

Design notes
------------
Files of this kind can be very large (RNA12 has ~1.6e7 genotypes, and
permissible sets average ~45-130 structures each, so a phenotype-grouped
file can have hundreds of millions of genotype references). To keep
memory at O(n_genotypes) rather than O(total entries), this script makes
two streaming passes over permissible_sets.txt:

  Pass 1: accumulate, per genotype, the sum of propensities of all scored
          phenotypes in its permissible set (a single float array of
          length n_genotypes -- a few hundred MB at most).
  Pass 2: re-read the file, divide each phenotype's propensity by the
          genotype's accumulated sum, and stream the result straight to
          disk (never held fully in memory).

If this is still too slow on the largest (most promiscuous) rules,
consider swapping the `line.split()` + `np.array(..., dtype=...)` parsing
for `np.fromstring(remainder, dtype=..., sep=" ")` (faster but deprecated
in recent numpy), or precomputing a binary/memory-mapped version of
permissible_sets.txt.
"""

import sys
import time
import argparse
import numpy as np


def count_lines(path):
    n = 0
    with open(path, "rb") as f:
        for _ in f:
            n += 1
    return n


def load_propensities(path):
    """Return dict: phenotype (str) -> propensity (float)."""
    props = {}
    with open(path) as f:
        for line in f:
            parts = line.split()
            if not parts:
                continue
            phenotype, prop = parts[0], float(parts[1])
            props[phenotype] = prop
    return props


def compute_occupancy(permissible_sets_path, genotypes_path,
                       propensities_path, output_path,
                       id_dtype=np.int32, val_dtype=np.float64,
                       progress_every=20):
    t0 = time.time()

    print("Counting genotypes...", file=sys.stderr)
    n_genotypes = count_lines(genotypes_path)
    print(f"  {n_genotypes:,} genotypes", file=sys.stderr)

    print("Loading propensities...", file=sys.stderr)
    props = load_propensities(propensities_path)
    print(f"  {len(props):,} phenotypes with known propensity",
          file=sys.stderr)

    # ---- Pass 1: accumulate per-genotype propensity sums -----------------
    print("Pass 1/2: accumulating per-genotype propensity sums...",
          file=sys.stderr)
    sums = np.zeros(n_genotypes, dtype=val_dtype)

    n_lines = 0
    n_skipped = 0
    n_entries = 0
    with open(permissible_sets_path) as f:
        for line in f:
            n_lines += 1
            parts = line.split()
            if not parts:
                continue
            phenotype = parts[0]
            prop = props.get(phenotype)
            if prop is None:
                n_skipped += 1
                continue
            ids = np.array(parts[1:], dtype=id_dtype)
            n_entries += ids.size
            # np.add.at is correct even if a gid appears more than once
            # on the same line (shouldn't happen, but stays safe).
            np.add.at(sums, ids, prop)

            if n_lines % progress_every == 0:
                print(f"  pass 1: {n_lines} phenotype lines, "
                      f"{n_entries:,} entries, "
                      f"{time.time()-t0:.0f}s elapsed", file=sys.stderr)

    print(f"Pass 1 done: {n_lines} lines ({n_skipped} skipped - no known "
          f"propensity), {n_entries:,} total entries, "
          f"{time.time()-t0:.0f}s elapsed", file=sys.stderr)

    n_zero = int(np.sum(sums == 0))
    if n_zero:
        print(f"  note: {n_zero:,} genotypes never appear with a scored "
              f"phenotype and will be absent from every output line.",
              file=sys.stderr)

    # ---- Pass 2: compute occupancy and stream to output -------------------
    print("Pass 2/2: computing occupancy and writing output...",
          file=sys.stderr)
    n_lines2 = 0
    with open(permissible_sets_path) as fin, open(output_path, "w") as fout:
        for line in fin:
            n_lines2 += 1
            parts = line.split()
            if not parts:
                continue
            phenotype = parts[0]
            prop = props.get(phenotype)
            if prop is None:
                continue
            ids = np.array(parts[1:], dtype=id_dtype)
            occ = prop / sums[ids]
            tokens = (f"{g}:{o:.6g}" for g, o in zip(ids.tolist(),
                                                       occ.tolist()))
            fout.write(phenotype + " " + " ".join(tokens) + "\n")

            if n_lines2 % progress_every == 0:
                print(f"  pass 2: {n_lines2} phenotype lines, "
                      f"{time.time()-t0:.0f}s elapsed", file=sys.stderr)

    print(f"Done. Wrote {output_path}. Total time: "
          f"{time.time()-t0:.0f}s", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--permissible-sets", default="permissible_sets.txt")
    ap.add_argument("--genotypes", default="genotypes.txt")
    ap.add_argument("--propensities", default="mfe_propensities.csv")
    ap.add_argument("--output", default="phenotype_occupancy.txt")
    ap.add_argument("--id-dtype", default="int32",
                     choices=["int32", "int64"],
                     help="int32 is enough for up to ~2.1e9 genotype ids "
                          "(RNA12 has ~1.6e7) and roughly halves memory "
                          "use vs int64.")
    args = ap.parse_args()

    id_dtype = np.int32 if args.id_dtype == "int32" else np.int64

    compute_occupancy(args.permissible_sets, args.genotypes,
                       args.propensities, args.output, id_dtype=id_dtype)


if __name__ == "__main__":
    main()
