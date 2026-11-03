#!/usr/bin/env python3
"""
calibrate_diffusion.py

Fit the diffusion field's correlation strength (n_iterations, alpha) to
match the real spatial correlation of the true unfolded genotype set from
your ViennaRNA canonical GP map, so the synthetic diffusion-based
"unfolded" set (from diffusion_unfolded_field.py + threshold_field.py) is
comparably correlated, not just comparable in overall fraction.

Input file format
-------------------
A plain text file of unfolded genotype strings, one per line, all the
same length (12). The alphabet used in this file does NOT need to match
the "JKLM" alphabet used elsewhere -- it's auto-detected from whatever 4
distinct characters actually appear (e.g. "ACGU" for a real canonical
ViennaRNA export). This works because the calibration only compares
summary correlation statistics between the real and synthetic data, never
literal genotype identities, so the two alphabets never need to line up.

Calibration statistic
-----------------------
neighbor correlation: the Pearson correlation between a genotype's binary
unfolded/folded indicator and its neighbor's indicator (averaged over all
12 sequence positions for stability), computed once on the real data and
then, for a grid of (alpha, n_iterations) settings, on the synthetic
field thresholded to match the real data's overall unfolded fraction.
This is scale-invariant and directly comparable between the two
alphabets/datasets. The (alpha, n_iterations) combination whose synthetic
neighbor correlation is closest to the real one is reported as the best
fit.

For efficiency, diffusion is run once per alpha value up to the largest
requested iteration count, checkpointing the correlation at each
intermediate iteration along the way -- rather than restarting the
diffusion from scratch for every (alpha, T) pair.

Usage
------
    python calibrate_diffusion.py true_unfolded.txt --max-iterations 10 \\
        --alphas 0.3,0.5,0.7,0.9 --seed 0
"""

import sys
import time
import argparse
import numpy as np

from diffusion_unfolded_field import (
    N_GENOTYPES, LENGTH, neighbor_sum, neighbor_correlation
)


def load_true_unfolded_indicator(path, length=LENGTH):
    """Read a file of unfolded genotype strings (one per line) and build
    a full-length binary indicator array (1=unfolded, 0=folded) over the
    entire genotype space, auto-detecting the alphabet from the
    characters that appear in the file.

    Returns: (indicator, alphabet, encode_fn)
        indicator: float64 array of length alphabet_size**length (4**12)
        alphabet: detected alphabet string, sorted, e.g. "ACGU"
        encode_fn: seq -> integer id under this file's own encoding
            (leftmost character = least-significant base-4 digit, same
            convention as seq_to_id in diffusion_unfolded_field.py)
    """
    with open(path) as f:
        genotypes = [line.strip() for line in f if line.strip()]
    if not genotypes:
        raise ValueError(f"no genotypes found in {path}")

    observed = set()
    for g in genotypes:
        if len(g) != length:
            raise ValueError(f"genotype {g!r} has length {len(g)}, expected {length}")
        observed.update(g)
    alphabet = "".join(sorted(observed))
    if len(alphabet) != 4:
        raise ValueError(
            f"expected a 4-letter alphabet in {path}, found {len(alphabet)}: {alphabet!r}")

    letter_to_digit = {c: i for i, c in enumerate(alphabet)}

    def encode(seq):
        gid = 0
        for p, c in enumerate(seq):
            gid += letter_to_digit[c] * (4 ** p)
        return gid

    n_total = 4 ** length
    indicator = np.zeros(n_total, dtype=np.float64)
    for g in genotypes:
        indicator[encode(g)] = 1.0

    return indicator, alphabet, encode


def avg_neighbor_correlation(indicator, n_positions=LENGTH):
    """Average neighbor_correlation over all sequence positions, for a
    more stable single summary statistic than any one position alone."""
    return float(np.mean([
        neighbor_correlation(indicator, position=p) for p in range(n_positions)
    ]))


def calibrate(true_path, max_iterations=10, alphas=(0.3, 0.5, 0.7, 0.9),
              seed=0, checkpoint_every=1, corr_positions=4, progress=True):
    """Sweep (alpha, n_iterations) and report how closely each setting's
    synthetic neighbor correlation matches the real data's. `corr_positions`
    controls how many of the 12 positions are averaged for the (repeated,
    per-checkpoint) synthetic correlation estimate -- fewer positions is
    faster but noisier; the real data's own correlation is always
    averaged over all 12 positions since it's only computed once.

    Returns: (results, real_corr, target_fraction)
        results: list of dicts sorted by closeness to the real
            correlation, each with keys alpha, n_iterations,
            synthetic_corr, target_corr, abs_diff, achieved_fraction
    """
    real_indicator, alphabet, _ = load_true_unfolded_indicator(true_path)
    target_fraction = float(real_indicator.mean())
    real_corr = avg_neighbor_correlation(real_indicator, n_positions=LENGTH)

    if progress:
        print(f"Real data: alphabet={alphabet!r}, "
              f"unfolded fraction={target_fraction:.4f}, "
              f"avg neighbor correlation={real_corr:.4f}", file=sys.stderr)

    results = []
    t0 = time.time()
    for alpha in alphas:
        rng = np.random.default_rng(seed)
        Z = rng.standard_normal(N_GENOTYPES)
        for T in range(1, max_iterations + 1):
            nsum = neighbor_sum(Z)
            Z = (1 - alpha) * Z + alpha * (nsum / (3 * LENGTH))
            if T % checkpoint_every != 0:
                continue
            threshold = np.quantile(Z, 1 - target_fraction)
            synthetic_indicator = (Z > threshold).astype(np.float64)
            achieved_fraction = float(synthetic_indicator.mean())
            synth_corr = avg_neighbor_correlation(
                synthetic_indicator, n_positions=corr_positions)
            results.append({
                "alpha": alpha,
                "n_iterations": T,
                "synthetic_corr": synth_corr,
                "target_corr": real_corr,
                "abs_diff": abs(synth_corr - real_corr),
                "achieved_fraction": achieved_fraction,
            })
            if progress:
                print(f"  alpha={alpha} T={T}: synthetic corr={synth_corr:.4f} "
                      f"(target {real_corr:.4f}), fraction={achieved_fraction:.4f} "
                      f"[{time.time()-t0:.0f}s elapsed]", file=sys.stderr)

    results.sort(key=lambda r: r["abs_diff"])
    return results, real_corr, target_fraction


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("true_unfolded_file",
                     help="Plain text file of true unfolded genotype "
                          "strings, one per line (any 4-letter alphabet).")
    ap.add_argument("--max-iterations", type=int, default=10)
    ap.add_argument("--alphas", default="0.3,0.5,0.7,0.9",
                     help="Comma-separated list of alpha values to sweep.")
    ap.add_argument("--checkpoint-every", type=int, default=1)
    ap.add_argument("--corr-positions", type=int, default=4,
                     help="Number of sequence positions to average for "
                          "the (repeated) synthetic correlation estimate "
                          "at each checkpoint -- lower is faster.")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--top", type=int, default=5,
                     help="How many best-matching settings to print.")
    args = ap.parse_args()

    alphas = tuple(float(a) for a in args.alphas.split(","))

    results, real_corr, target_fraction = calibrate(
        args.true_unfolded_file, max_iterations=args.max_iterations,
        alphas=alphas, seed=args.seed, checkpoint_every=args.checkpoint_every,
        corr_positions=args.corr_positions)

    print(f"\nTarget: unfolded fraction={target_fraction:.4f}, "
          f"neighbor correlation={real_corr:.4f}\n")
    print(f"Best-matching settings (closest synthetic correlation first):")
    for r in results[:args.top]:
        print(f"  alpha={r['alpha']}, n_iterations={r['n_iterations']}: "
              f"synthetic corr={r['synthetic_corr']:.4f} "
              f"(diff {r['abs_diff']:.4f}), "
              f"achieved fraction={r['achieved_fraction']:.4f}")


if __name__ == "__main__":
    main()
