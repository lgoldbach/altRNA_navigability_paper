#!/usr/bin/env python3
"""
select_unfolded_components.py

Randomly select whole neutral components (networks) -- never split one --
so that the total number of genotypes they cover is close to a target
count N, with a margin of error (an exact match is generally impossible,
since you can only add or omit a component as a whole unit; the closest
achievable total depends on the size of whichever component happens to be
the last one considered).

Input: nc_to_gt, a dict {component_id: [genotype_string, ...]}.

Two ways to pick components, controlled by `mode`:
  "size_weighted" (default) -- components are drawn without replacement
      with probability proportional to their size, so bigger components
      tend to be picked earlier. This is the direct network-respecting
      analogue of "repeatedly pick a uniformly random genotype and mark
      its whole component" -- i.e. it mirrors the original per-genotype
      random-assignment procedure as closely as possible while still
      keeping components intact.
  "uniform" -- every component is equally likely to be picked next,
      regardless of size.

Algorithm: shuffle the components (weighted or uniform), accumulate them
in that order until the running total first reaches or exceeds the
target, then check whether stopping one component earlier (i.e. omitting
the one that pushed the total over) would have landed closer to the
target -- take whichever of the two is closer. If the result isn't within
`margin` of the target, reshuffle and try again (up to max_attempts),
keeping track of the best attempt seen so we always return *something*
usable even if no attempt lands inside the margin.
"""

import math
import random


def _weighted_shuffle(ids, sizes, rng):
    """Random ordering without replacement, with selection probability at
    each step proportional to remaining size (Efraimidis-Spirakis order
    sampling). Done in log-space so it stays numerically stable even for
    very large components (a naive u**(1/w) key collapses to ~1.0 for
    large w and loses all discriminating precision)."""
    keyed = []
    for cid in ids:
        u = rng.random()
        w = sizes[cid]
        log_key = float("-inf") if w <= 0 else math.log(u) / w
        keyed.append((log_key, cid))
    keyed.sort(reverse=True)
    return [cid for _, cid in keyed]


def _uniform_shuffle(ids, rng):
    order = list(ids)
    rng.shuffle(order)
    return order


def select_unfolded_components(nc_to_gt, target_n, margin=None,
                                mode="size_weighted", max_attempts=1000,
                                seed=None):
    """
    nc_to_gt: dict {component_id: [genotype_string, ...]}
    target_n: desired total number of genotypes across selected components
    margin: acceptable absolute deviation |achieved_n - target_n|.
        Default: max(1% of target_n, size of the largest component) --
        the latter floor matters because you fundamentally cannot land
        closer than "the size of whichever single component tips the
        balance" if that component is large, so a margin smaller than the
        biggest component may simply be unreachable.
    mode: "size_weighted" (default) or "uniform" -- see module docstring.
    max_attempts: number of random orderings to try before giving up and
        returning the best (closest-to-target) one found.
    seed: optional int for reproducibility.

    Returns a dict:
        selected_ids   -- list of component ids included
        genotypes      -- set of genotype strings covered by those ids
        achieved_n     -- len(genotypes)
        target_n       -- the requested target (echoed back)
        deviation       -- achieved_n - target_n (signed)
        within_margin  -- bool, whether |deviation| <= margin
        margin         -- the margin actually used
        attempts_used  -- how many random orderings were tried
    """
    rng = random.Random(seed)
    ids = list(nc_to_gt.keys())
    sizes = {cid: len(nc_to_gt[cid]) for cid in ids}
    total = sum(sizes.values())

    if target_n < 0:
        raise ValueError("target_n must be >= 0")
    if target_n > total:
        raise ValueError(
            f"target_n ({target_n}) exceeds total genotypes available "
            f"({total}) across all {len(ids)} components")

    if margin is None:
        margin = max(int(round(0.01 * target_n)), max(sizes.values(), default=0))

    if mode not in ("size_weighted", "uniform"):
        raise ValueError(f"unknown mode {mode!r}, use 'size_weighted' or 'uniform'")

    if target_n == 0:
        return {
            "selected_ids": [],
            "genotypes": set(),
            "achieved_n": 0,
            "target_n": 0,
            "deviation": 0,
            "within_margin": True,
            "margin": margin,
            "attempts_used": 0,
        }

    best_selected = []
    best_total = 0
    best_deviation = -target_n  # deviation of the empty selection

    attempts_used = 0
    for attempt in range(1, max_attempts + 1):
        attempts_used = attempt
        order = (_weighted_shuffle(ids, sizes, rng) if mode == "size_weighted"
                  else _uniform_shuffle(ids, rng))

        selected = []
        running_total = 0
        for cid in order:
            selected.append(cid)
            running_total += sizes[cid]
            if running_total >= target_n:
                break

        # stopping one component earlier might land closer to target
        if len(selected) > 1:
            without_last = running_total - sizes[selected[-1]]
            if abs(without_last - target_n) < abs(running_total - target_n):
                selected = selected[:-1]
                running_total = without_last

        deviation = running_total - target_n

        if abs(deviation) < abs(best_deviation):
            best_selected, best_total, best_deviation = selected, running_total, deviation

        if abs(deviation) <= margin:
            best_selected, best_total, best_deviation = selected, running_total, deviation
            break

    genotypes = set()
    for cid in best_selected:
        genotypes.update(nc_to_gt[cid])

    return {
        "selected_ids": best_selected,
        "genotypes": genotypes,
        "achieved_n": best_total,
        "target_n": target_n,
        "deviation": best_deviation,
        "within_margin": abs(best_deviation) <= margin,
        "margin": margin,
        "attempts_used": attempts_used,
    }


if __name__ == "__main__":
    # small demo
    rng = random.Random(0)
    nc_to_gt = {
        i: [f"g{i}_{j}" for j in range(rng.choice([1, 2, 5, 20, 100, 1000, 50000]))]
        for i in range(200)
    }
    total = sum(len(v) for v in nc_to_gt.values())
    print(f"total genotypes across {len(nc_to_gt)} components: {total}")

    target = int(0.6 * total)
    result = select_unfolded_components(nc_to_gt, target, seed=1)
    print(f"target={target}, achieved={result['achieved_n']}, "
          f"deviation={result['deviation']}, margin={result['margin']}, "
          f"within_margin={result['within_margin']}, "
          f"attempts_used={result['attempts_used']}, "
          f"n_components_selected={len(result['selected_ids'])}")
