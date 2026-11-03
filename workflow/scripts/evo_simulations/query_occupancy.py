#!/usr/bin/env python3
"""
query_occupancy.py

Look up the propensity-weighted occupancy of one target ("functional")
phenotype for every genotype, from the phenotype_occupancy.txt file
produced by compute_occupancy.py, and optionally turn it into a dense
per-genotype fitness/occupancy vector for use directly in landscape /
adaptive-walk code.

Genotypes not listed for the target phenotype implicitly have occupancy 0
(the phenotype is not in their permissible set).

Library usage
-------------
    from query_occupancy import load_occupancy, occupancy_vector

    occ = load_occupancy("phenotype_occupancy.txt", "(((...)))")
    # occ: dict {genotype_id: occupancy}, only nonzero entries

    fitness = occupancy_vector(occ, n_genotypes=16777216)
    # dense numpy array, length n_genotypes, 0 everywhere except occ's keys

Command-line usage
-------------------
    python query_occupancy.py phenotype_occupancy.txt "(((...)))" \\
        --n-genotypes 16777216 --out fitness_target.npy
"""

import sys
import argparse
import numpy as np


def load_occupancy(path, target_phenotype):
    """Scan phenotype_occupancy.txt for the target phenotype's line and
    return {genotype_id: occupancy}. The file has one line per scored
    phenotype (typically tens to a few hundred lines), so a linear scan
    is fast regardless of how many genotypes are listed per line."""
    with open(path) as f:
        for line in f:
            parts = line.split()
            if not parts or parts[0] != target_phenotype:
                continue
            occ = {}
            for tok in parts[1:]:
                gid_str, occ_str = tok.split(":")
                occ[int(gid_str)] = float(occ_str)
            return occ
    raise KeyError(
        f"Phenotype {target_phenotype!r} not found in {path} (it may lack "
        f"a known propensity, or never occur in any genotype's "
        f"permissible set)."
    )


def occupancy_vector(occ, n_genotypes, dtype=np.float64):
    """Turn a {genotype_id: occupancy} dict into a dense length-n_genotypes
    numpy array (0 for every genotype not in occ)."""
    vec = np.zeros(n_genotypes, dtype=dtype)
    if occ:
        ids = np.fromiter(occ.keys(), dtype=np.int64, count=len(occ))
        vals = np.fromiter(occ.values(), dtype=dtype, count=len(occ))
        vec[ids] = vals
    return vec


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("occupancy_file")
    ap.add_argument("target_phenotype")
    ap.add_argument("--n-genotypes", type=int, required=True,
                     help="Total number of genotypes (line count of "
                          "genotypes.txt), needed to size the dense "
                          "output vector.")
    ap.add_argument("--out", default=None,
                     help="Optional .npy path to save the dense "
                          "fitness/occupancy vector.")
    args = ap.parse_args()

    occ = load_occupancy(args.occupancy_file, args.target_phenotype)
    print(f"{len(occ):,} genotypes have nonzero occupancy for "
          f"{args.target_phenotype!r}", file=sys.stderr)

    if args.out:
        vec = occupancy_vector(occ, args.n_genotypes)
        np.save(args.out, vec)
        print(f"Saved dense vector ({vec.size:,} entries) to {args.out}",
              file=sys.stderr)


if __name__ == "__main__":
    main()
