#!/usr/bin/env python3
"""
similarity_fitness_landscape.py

Generate correlated ("single-peak") fitness landscapes in which one
randomly chosen phenotype is the fittest, and every other phenotype's
fitness decays with its structural (base-pair) distance to that target.
This replaces the paper's i.i.d.-Uniform(0,1) fitness assignment with a
minimal, tunable form of correlated landscape, as a robustness check
against the independence assumption.

Distance metric
----------------
Two phenotypes are compared by base-pair distance: the size of the
symmetric difference between their base-pair sets. This only needs the
dot-bracket string (no ViennaRNA / RNAdistance call needed), works for
non-canonical structures the same as canonical ones, and is well defined
because all phenotypes here are secondary structures of RNAs of the same
length. It is the standard structure-distance measure in this literature
(what RNAdistance's -Dp mode computes).

Fitness assignment
-------------------
Given target T and distance d(S, T) for phenotype S, a deterministic decay
component:
    "exponential": f_decay(S) = exp(-d(S,T) / corr_length)
    "linear":      f_decay(S) = max(0, 1 - d(S,T) / max_d)

Both are bounded in [0, 1] with f_decay(T) = 1. corr_length controls how
quickly f_decay falls off with distance: *small* corr_length -> a sharp
spike (only structures nearly identical to T get appreciable fitness,
everything else ~0); *large* corr_length -> f_decay flattens out toward 1
for essentially everything (not toward 0 -- the exponent -d/corr_length
shrinks to 0 as corr_length grows, so exp(...) -> 1).

Importantly, f_decay alone is purely deterministic given the target and
has no regime that recovers the paper's original i.i.d.-Uniform(0,1)
landscape: the small-corr_length limit is a degenerate single-genotype
spike, and the large-corr_length limit is a flat, uninformative landscape
close to 1 everywhere. Neither is "uncorrelated." To get an actual
correlated <-> uncorrelated sweep, fitness is a mix of f_decay and an
independent random draw per phenotype:

    fitness(S) = mix * f_decay(S) + (1 - mix) * U(S),   U(S) ~ Uniform(0,1) i.i.d.

This is the standard "Rough Mount Fuji" construction (Aita et al. 2000;
Neidhart, Szendro & Krug 2014): a deterministic, distance-based component
plus i.i.d. noise, combined with weight `mix`.
    mix = 0 -> fitness(S) = U(S) for every phenotype, independently drawn
               -- this is *exactly* the paper's original model (not just
               close to it), since U(S) is i.i.d. Uniform(0,1) per
               phenotype with no reference to the target at all.
    mix = 1 -> purely the deterministic single-peak landscape described
               above (no randomness).
    0 < mix < 1 -> a genuine blend, and the parameter to sweep for the
               robustness analysis: how much correlation can be
               introduced before the base-pairing-rule ranking changes.

Target selection is uniform random over phenotypes (not weighted by
phenotype frequency), matching how the paper's original model gives every
phenotype an equal chance of being the global peak.

Usage as a library
-------------------
    from similarity_fitness_landscape import (
        load_phenotypes, generate_landscapes, write_landscapes_separate)

    phenotypes = load_phenotypes("mfe_propensities.csv")  # or any file
                                                            # with phenotype
                                                            # as first token
    landscapes = generate_landscapes(phenotypes, n_landscapes=200,
                                      decay="exponential", corr_length=3.0,
                                      mix=1.0, seed=0)
    target, fitness = landscapes[0]
    # fitness: dict {phenotype: fitness value in [0, 1]}

    # one .txt file per landscape (phenotype fitness, one per line),
    # plus a landscape_targets.txt manifest:
    write_landscapes_separate(landscapes, "similarity_landscapes")

Command-line usage
-------------------
    python similarity_fitness_landscape.py mfe_propensities.csv \\
        --n-landscapes 200 --decay exponential --corr-length 3.0 --mix 1.0 \\
        --seed 0 --output-dir similarity_landscapes
    # -> similarity_landscapes/landscape_0000.txt, landscape_0001.txt, ...
    #    similarity_landscapes/landscape_targets.txt
"""

import sys
import random
import argparse
import numpy as np


def parse_pairs(dotbracket):
    """Return the set of base pairs (i, j), i < j, in a pseudoknot-free
    dot-bracket string, via a simple stack (matches the paper's folding
    constraints: nested pairs only)."""
    stack = []
    pairs = set()
    for i, c in enumerate(dotbracket):
        if c == "(":
            stack.append(i)
        elif c == ")":
            j = stack.pop()
            pairs.add((j, i))
    return pairs


def base_pair_distance(pairs1, pairs2):
    """Symmetric-difference size between two base-pair sets."""
    return len(pairs1 ^ pairs2)


def load_phenotypes(path):
    """Read the unique phenotype (first whitespace-separated token) from
    every line of a file -- works directly on mfe_propensities.csv, or on
    permissible_sets.txt, or on a plain one-phenotype-per-line file."""
    phenotypes = []
    seen = set()
    with open(path) as f:
        for line in f:
            parts = line.split()
            if not parts:
                continue
            p = parts[0]
            if p not in seen:
                seen.add(p)
                phenotypes.append(p)
    return phenotypes


def build_landscape(phenotypes, pairs_cache, target=None, decay="exponential",
                     corr_length=3.0, mix=1.0, rng=None):
    """Build one similarity-based fitness landscape.

    phenotypes: list of unique phenotype strings
    pairs_cache: dict phenotype -> base-pair set (precomputed once, reused
        across many landscapes for speed)
    target: force a specific target phenotype instead of drawing randomly
    mix: weight in [0, 1] on the deterministic distance-decay component
        vs. an i.i.d. Uniform(0,1) random component drawn independently
        per phenotype:
            fitness(S) = mix * f_decay(S) + (1 - mix) * U(S)
        mix=1.0 (default) -> purely deterministic single-peak landscape.
        mix=0.0 -> exactly the paper's original i.i.d. Uniform(0,1)
                   landscape (f_decay/target play no role at all).
    rng: a random.Random instance (for reproducibility across repeated
        calls); defaults to the global random module
    Returns: (target_phenotype, {phenotype: fitness})
    """
    rng = rng or random
    if target is None:
        target = rng.choice(phenotypes)

    target_pairs = pairs_cache[target]
    dist = {p: base_pair_distance(pairs_cache[p], target_pairs)
            for p in phenotypes}

    if decay == "exponential":
        f_decay = {p: float(np.exp(-d / corr_length)) for p, d in dist.items()}
    elif decay == "linear":
        d_max = max(dist.values()) or 1
        f_decay = {p: max(0.0, 1.0 - d / d_max) for p, d in dist.items()}
    else:
        raise ValueError(f"unknown decay {decay!r}, use 'exponential' or 'linear'")

    if mix >= 1.0:
        fitness = f_decay
    else:
        u = {p: rng.random() for p in phenotypes}
        if mix <= 0.0:
            fitness = u
        else:
            fitness = {p: mix * f_decay[p] + (1 - mix) * u[p] for p in phenotypes}

    return target, fitness


def generate_landscapes(phenotypes, n_landscapes=200, decay="exponential",
                         corr_length=3.0, mix=1.0, seed=0):
    """Generate n_landscapes replicate similarity-based landscapes, each
    with an independently, uniformly drawn random target phenotype --
    the direct analogue of the paper's 200-random-uniform-landscapes
    setup. `mix` is the same correlated<->uncorrelated blend weight as in
    build_landscape (mix=0 reproduces the paper's original i.i.d. model
    exactly). Returns a list of (target, {phenotype: fitness}) tuples."""
    rng = random.Random(seed)
    pairs_cache = {p: parse_pairs(p) for p in phenotypes}
    landscapes = []
    for _ in range(n_landscapes):
        target, fitness = build_landscape(phenotypes, pairs_cache,
                                           decay=decay,
                                           corr_length=corr_length,
                                           mix=mix, rng=rng)
        landscapes.append((target, fitness))
    return landscapes


def write_landscapes(landscapes, output_path):
    """Write landscapes to a single text file, one line per landscape:
    <target_phenotype> <phenotype_1>:<fitness_1> <phenotype_2>:<fitness_2> ...
    (same 'token:value' convention as phenotype_occupancy.txt)."""
    with open(output_path, "w") as f:
        for target, fitness in landscapes:
            tokens = (f"{p}:{v:.6g}" for p, v in fitness.items())
            f.write(target + " " + " ".join(tokens) + "\n")


def load_landscapes(path):
    """Read landscapes back from the file written by write_landscapes."""
    landscapes = []
    with open(path) as f:
        for line in f:
            parts = line.split()
            target = parts[0]
            fitness = {}
            for tok in parts[1:]:
                p, v = tok.rsplit(":", 1)
                fitness[p] = float(v)
            landscapes.append((target, fitness))
    return landscapes


def write_landscapes_separate(landscapes, output_dir, prefix="landscape",
                               digits=None):
    """Write one .txt file per landscape into output_dir, each containing
    one phenotype per line followed by its fitness, space-separated
    (same layout as mfe_propensities.csv):

        (((...))) 0.842
        ........... 0.113
        ((.....)) 0.5

    Files are named <prefix>_<index>.txt, zero-padded (e.g.
    landscape_0000.txt, landscape_0001.txt, ...).

    Since fitness(target) isn't necessarily 1 once mix < 1 (the random
    component can push it down), the target phenotype can't always be
    read back off the landscape file by looking for the max value alone.
    A separate manifest <prefix>_targets.txt is written alongside the
    per-landscape files, one line per landscape: "<index> <target>" --
    this keeps every landscape file in the plain two-column format you
    asked for, with target metadata kept out of it.

    Returns the list of landscape file paths written, in order.
    """
    import os

    os.makedirs(output_dir, exist_ok=True)
    n = len(landscapes)
    digits = digits or max(4, len(str(max(n - 1, 0))))

    paths = []
    targets_path = os.path.join(output_dir, f"{prefix}_targets.txt")
    with open(targets_path, "w") as tf:
        for i, (target, fitness) in enumerate(landscapes):
            fname = f"{prefix}_{i:0{digits}d}.txt"
            fpath = os.path.join(output_dir, fname)
            with open(fpath, "w") as f:
                for p, v in fitness.items():
                    f.write(f"{p} {v:.6g}\n")
            tf.write(f"{i} {target}\n")
            paths.append(fpath)

    return paths


def load_landscapes_separate(output_dir, prefix="landscape"):
    """Read landscapes back from the files written by
    write_landscapes_separate. Returns a list of (target, {phenotype:
    fitness}) tuples, in the same order as the targets manifest."""
    import os

    targets_path = os.path.join(output_dir, f"{prefix}_targets.txt")
    entries = []
    with open(targets_path) as tf:
        for line in tf:
            parts = line.split()
            if not parts:
                continue
            idx, target = parts[0], parts[1]
            entries.append((int(idx), target))

    landscapes = []
    for idx, target in entries:
        # zero-padding width isn't stored explicitly, so just try widths
        # until the file is found
        fpath = None
        for width in range(1, 9):
            candidate = os.path.join(output_dir, f"{prefix}_{idx:0{width}d}.txt")
            if os.path.exists(candidate):
                fpath = candidate
                break
        if fpath is None:
            raise FileNotFoundError(
                f"no landscape file found for index {idx} in {output_dir}")

        fitness = {}
        with open(fpath) as f:
            for line in f:
                parts = line.split()
                if not parts:
                    continue
                fitness[parts[0]] = float(parts[1])
        landscapes.append((target, fitness))

    return landscapes


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("phenotype_file",
                     help="Any file with the phenotype as the first "
                          "whitespace-separated token per line, e.g. "
                          "mfe_propensities.csv or permissible_sets.txt")
    ap.add_argument("--n-landscapes", type=int, default=200)
    ap.add_argument("--decay", choices=["exponential", "linear"],
                     default="exponential")
    ap.add_argument("--corr-length", type=float, default=3.0,
                     help="Only used for --decay exponential. Smaller = "
                          "sharper single spike near the target; larger = "
                          "flatter landscape close to 1 everywhere. "
                          "Neither extreme is 'uncorrelated' -- use --mix "
                          "for that.")
    ap.add_argument("--mix", type=float, default=1.0,
                     help="Weight in [0,1] on the deterministic decay vs. "
                          "an i.i.d. Uniform(0,1) random component per "
                          "phenotype. mix=1 -> pure similarity landscape "
                          "(default). mix=0 -> exactly the paper's "
                          "original i.i.d. random landscape. Sweep this "
                          "for the actual correlated<->uncorrelated "
                          "robustness check.")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--output-dir", default="similarity_landscapes",
                     help="Directory to write one .txt file per landscape "
                          "into (phenotype fitness, one per line, plus a "
                          "<prefix>_targets.txt manifest). This is the "
                          "default output mode.")
    ap.add_argument("--prefix", default="landscape",
                     help="Filename prefix for the per-landscape files.")
    ap.add_argument("--single-file", default=None,
                     help="If given, ALSO write all landscapes into one "
                          "combined file at this path (phenotype:fitness "
                          "tokens per line), in addition to the "
                          "per-landscape files in --output-dir.")
    args = ap.parse_args()

    phenotypes = load_phenotypes(args.phenotype_file)
    print(f"{len(phenotypes)} unique phenotypes loaded from "
          f"{args.phenotype_file}", file=sys.stderr)

    landscapes = generate_landscapes(
        phenotypes, n_landscapes=args.n_landscapes, decay=args.decay,
        corr_length=args.corr_length, mix=args.mix, seed=args.seed)

    paths = write_landscapes_separate(landscapes, args.output_dir,
                                       prefix=args.prefix)
    print(f"Wrote {len(paths)} landscape files to {args.output_dir}/ "
          f"(plus {args.prefix}_targets.txt)", file=sys.stderr)

    if args.single_file:
        write_landscapes(landscapes, args.single_file)
        print(f"Also wrote combined file to {args.single_file}",
              file=sys.stderr)


if __name__ == "__main__":
    main()
