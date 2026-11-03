#!/usr/bin/env python3
"""
diffusion_unfolded_field.py

Generate a spatially correlated random field over the RNA12 genotype
space (length 12, alphabet size 4: J, K, L, M) via iterative graph
diffusion (smoothing) -- a correlated alternative to i.i.d. random
genotype selection, reproducing realistic local correlation on the
mutational (Hamming) graph.

Genotype space and neighbor structure
---------------------------------------
Genotypes are length-12 strings over "JKLM". Each has exactly 36
single-mutation neighbors (12 positions x 3 alternative letters). Rather
than build an explicit neighbor list (36 x 4^12 ~ 6x10^8 edges), this
exploits the structure of the genotype space directly: it's a
12-dimensional grid of side 4 (a generalized hypercube / Hamming graph),
so summing a field over "all genotypes differing from g at position p"
is just a sum along one axis of a 12-dimensional array -- computable for
every genotype at once via a single numpy reduction. Repeating over the
12 positions and adding up gives, for every genotype simultaneously, the
sum of its 36 neighbors' values, with no per-genotype neighbor
enumeration at all.

Encoding: genotype string s (s[0] leftmost .. s[11] rightmost), digit(c)
= "JKLM".index(c). id(s) = sum_p digit(s[p]) * 4**p, i.e. s[0] is the
least-significant base-4 digit. See seq_to_id / id_to_seq.

Algorithm
----------
1. Z(g) ~ i.i.d. N(0,1) for all 4^12 genotypes.
2. Repeat T times: Z(g) <- (1-alpha)*Z(g) + alpha*mean(Z over g's 36
   neighbors).
3. Result: a smooth, correlated field -- genotypes close together in
   mutation space end up with similar values.
4. Write every genotype and its final value to a file. Use
   threshold_field.py to pick a cutoff for a target "unfolded" fraction.

A note on tuning T and alpha: each genotype has 36 neighbors (a very
densely connected graph), so this mixes toward the global mean much
faster than 1D/2D diffusion intuition suggests -- empirically, the raw
field's standard deviation roughly halves per iteration at alpha=0.5, and
by T~15 the field is nearly constant everywhere. This does NOT hurt
threshold_field.py's usability (it selects by quantile, which only cares
about relative ordering, not absolute scale), but it does mean T should
be a handful of iterations (think low single digits), not dozens.
neighbor_correlation() below computes the exact correlation between
genotypes and their neighbors for a given field, scale-invariant and
directly comparable across settings -- use it to calibrate T/alpha
against whatever real neighbor-correlation statistic you measure from
your true (ViennaRNA-derived) canonical unfolded set, per the earlier
discussion, rather than guessing. For reference, at alpha=0.5: T=1 gives
neighbor corr ~0.06, T=4 ~0.22, T=8 ~0.40; at alpha=0.8: T=1 ~0.17, T=4
~0.43, T=8 ~0.64.

Usage
------
    python diffusion_unfolded_field.py --iterations 4 --alpha 0.5 \\
        --seed 0 --output field.txt
"""

import sys
import time
import argparse
import itertools
import numpy as np

ALPHABET = "JKLM"
LENGTH = 12
N_GENOTYPES = len(ALPHABET) ** LENGTH
_LETTER_TO_DIGIT = {c: i for i, c in enumerate(ALPHABET)}


def seq_to_id(seq):
    """Genotype string -> integer id (0 .. 4^12-1). seq[0] is the
    least-significant base-4 digit."""
    gid = 0
    for p, c in enumerate(seq):
        gid += _LETTER_TO_DIGIT[c] * (4 ** p)
    return gid


def id_to_seq(gid):
    """Integer id -> genotype string."""
    digits = []
    for _ in range(LENGTH):
        digits.append(gid % 4)
        gid //= 4
    return "".join(ALPHABET[d] for d in digits)


def _axis_for_position(p):
    """Axis of the (4,)*12 C-order reshape of the flat genotype array
    that corresponds to sequence position p. Derivation: for a flat
    array indexed by id = sum_p digit_p * 4**p reshaped with C order (row
    major, last axis fastest), reshaped[a_0,...,a_11] maps to flat index
    a_0*4^11 + a_1*4^10 + ... + a_11*4^0. Matching coefficients of 4^m to
    id's digit_m term gives digit_m <-> axis index (11-m)."""
    return LENGTH - 1 - p


def neighbor_sum(Z_flat):
    """Given Z as a flat length-4^12 array indexed by genotype id, return
    a flat array of the same length where entry g is the sum of Z over
    g's 36 single-mutation neighbors."""
    shape = (4,) * LENGTH
    Z = Z_flat.reshape(shape)
    total = np.zeros(shape, dtype=Z.dtype)
    for p in range(LENGTH):
        axis = _axis_for_position(p)
        axis_sum = Z.sum(axis=axis, keepdims=True)  # sum over all 4 values at this position
        total += axis_sum - Z  # drop the self-term, leaving the 3 alternatives
    return total.reshape(-1)


def neighbor_correlation(Z_flat, position=0):
    """Exact Pearson correlation between Z(g) and Z(g') for genotypes g,
    g' that are neighbors differing only at the given sequence position
    (cyclically shifted by +1 letter) -- computed for the whole genotype
    space at once via np.roll along the corresponding axis. Scale
    invariant, so useful for comparing correlation strength across
    different alpha/n_iterations settings regardless of how much the raw
    field's variance has shrunk. By symmetry of the graph this doesn't
    depend much on which position is chosen; position=0 is a reasonable
    default single-position proxy for overall neighbor correlation."""
    shape = (4,) * LENGTH
    Z = Z_flat.reshape(shape)
    axis = _axis_for_position(position)
    Z_shifted = np.roll(Z, shift=-1, axis=axis)
    return float(np.corrcoef(Z.reshape(-1), Z_shifted.reshape(-1))[0, 1])


def diffuse(n_iterations=20, alpha=0.5, seed=0, dtype=np.float64,
            progress=True):
    """Run the iterative diffusion process; return the final flat field
    (length 4^12, indexed by genotype id -- see seq_to_id)."""
    rng = np.random.default_rng(seed)
    Z = rng.standard_normal(N_GENOTYPES).astype(dtype)

    t0 = time.time()
    for it in range(1, n_iterations + 1):
        nsum = neighbor_sum(Z)
        Z = (1 - alpha) * Z + alpha * (nsum / (3 * LENGTH))
        if progress:
            print(f"  diffusion iteration {it}/{n_iterations} "
                  f"({time.time()-t0:.1f}s elapsed)", file=sys.stderr)
    return Z


def _sequences_in_id_order():
    """Yield genotype strings in increasing id order (0, 1, 2, ...),
    without ever materializing the id explicitly. itertools.product
    varies its LAST argument fastest; since seq[0] is the
    least-significant digit (fastest-varying as id increases), each
    product tuple must be reversed to become a genotype string."""
    for combo in itertools.product(ALPHABET, repeat=LENGTH):
        yield "".join(reversed(combo))


def write_field_file(Z_flat, output_path, chunk_size=500_000, progress=True):
    """Write 'genotype value' (space-separated), one line per genotype,
    in increasing id order, to output_path."""
    t0 = time.time()
    written = 0
    with open(output_path, "w") as f:
        buf = []
        for seq, val in zip(_sequences_in_id_order(), Z_flat):
            buf.append(f"{seq} {val:.6g}\n")
            written += 1
            if len(buf) >= chunk_size:
                f.write("".join(buf))
                buf = []
                if progress:
                    print(f"  wrote {written:,}/{N_GENOTYPES:,} genotypes "
                          f"({time.time()-t0:.1f}s elapsed)", file=sys.stderr)
        if buf:
            f.write("".join(buf))
    if progress:
        print(f"Done writing {written:,} genotypes to {output_path} "
              f"({time.time()-t0:.1f}s total)", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--iterations", type=int, default=4,
                     help="Number of diffusion iterations. Each genotype "
                          "has 36 neighbors, so this graph mixes fast -- "
                          "keep this in the low single digits unless you "
                          "have measured that you need more correlation; "
                          "see neighbor_correlation() to calibrate.")
    ap.add_argument("--alpha", type=float, default=0.5,
                     help="Smoothing weight per iteration, in [0,1]. "
                          "0 = no smoothing (stays i.i.d.), 1 = fully "
                          "replaced by neighbor mean each step.")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--output", default="field.txt")
    args = ap.parse_args()

    print(f"Generating diffusion field over {N_GENOTYPES:,} genotypes "
          f"(L={LENGTH}, alphabet={ALPHABET})...", file=sys.stderr)
    Z = diffuse(n_iterations=args.iterations, alpha=args.alpha, seed=args.seed)
    write_field_file(Z, args.output)


if __name__ == "__main__":
    main()
