#!/usr/bin/env python

import pickle
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import distance
from rna_gpf.parsing import many_to_one_map_from_file_to_dict
from rna_gpf.base_pairing import BasePairing
from rna_gpf.utils import random_mutation_at_sites
from matplotlib.ticker import MaxNLocator


def compute_ph_diff_prob_of_pairs_by_hamming_dist(gp_map, sample_size, rng):
    hdist_to_ph_diff_prob = {}
    for ham_dist in range(1, 11):  # for each hamming distance from 1-10
        hdist_to_ph_diff_prob[ham_dist] = 0
        for i in range(sample_size):  # get N genotypes at that distance
            gt1 = rng.choice(gp_map.genotypes)
            sites = rng.choice(range(12), size=ham_dist, replace=False)  # find random unique sites to mutate
            gt2_folded = False
            for attempt in range(100):  # 100 attempts at finding folded gt2
                gt2 = random_mutation_at_sites(gt1, sites, gp_map.alphabet)  # mutate into second gt
                try:  # mapping gt2 does not raise error, it is folded 
                    ph2 = gp_map.map(gt2)
                    gt2_folded = True
                except KeyError:  
                    continue  # if gt2 unfolded we try again

                if gt2_folded:  # if gt2 is folded we can continue
                    if gp_map.map(gt1) != ph2:  # check if phenotypes differ
                        hdist_to_ph_diff_prob[ham_dist] += 1
                    break  # done with attempts
        
        hdist_to_ph_diff_prob[ham_dist] /= sample_size
    
    return hdist_to_ph_diff_prob

red_map = {1: {'J': 'J', 'K': 'R', 'L': 'R', 'M': 'M'},
               2: {'J': 'R', 'K': 'R', 'L': 'L', 'M': 'M'},
               6: {'M': 'M', 'K': 'R', 'L': 'R', 'J': 'R'},
               7: {'J': 'X', 'K': 'X', 'L': 'Y', 'M': 'Y'},  # X and why are two different redundants
               9: {'J': 'R', 'K': 'R', 'L': 'L', 'M': 'M'}}

def identical_sites(gt1, gt2, l):
    d = l-round(distance.hamming(list(gt1), list(gt2)) * l)
    return d

def map_redund_sites_to_same_letter(gts, bp_rule):
    new_gts = []
    for gt in gts:
        new_gt = ""
        for l in gt:
            l_ = red_map[bp_rule][l]
            new_gt += l_
        new_gts.append(new_gt)
    return new_gts

def get_n_non_neutral_neighbors(gts, gp_map, n, unique_n, rng):
    non_neut_neighs = []
    for i in range(n):
        nc_gt = rng.choice(gts, size=1)[0]
        neighs = get_non_neutral_neighbors(gp_map, nc_gt)
        if neighs:
            non_neut_neigh = rng.choice(neighs, size=1)[0]
            non_neut_neighs.append(non_neut_neigh)
    non_neut_neighs = np.unique(non_neut_neighs)[:unique_n]
    print(f"Only found: {len(non_neut_neighs)} non-neutral neighbors")
    return non_neut_neighs

def 25.


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gp_map", help="GP map pickle files ", 
                        required=True, nargs='+')
    parser.add_argument("-n", "--nc_to_gt", help="neutral network to genotype map (.txt) ", 
                        required=True, nargs='+')
    parser.add_argument("-o", "--output", help="Output file name "
                        "(should end in .pdf)", required=True)

    args = parser.parse_args()
    fig, (ax1, ax2) = plt.subplots(ncols=2, nrows=1, figsize=(10, 5))
    
    ax1.set_box_aspect(1)
    ax2.set_box_aspect(1)

    id_sites_means = []
    deg_sites_means = []

    xs = []
    ys = []

    seed = np.random.randint(low=1, high=10000)
    # seed=
    rng = np.random.default_rng(seed)
    print(f"RANDOM SEED: {seed}")

    gave_label = False
    for bp, (gp_map_f, nc_to_gt_f) in enumerate(zip(args.gp_map, args.nc_to_gt), start=1):
        print(f"GP MAP {bp}")
        gp_map = pickle.load(open(gp_map_f, "rb"))

        hdist_to_ph_diff = compute_ph_diff_prob_of_pairs_by_hamming_dist(gp_map, sample_size=100, rng=rng)
        x, y = list(zip(*hdist_to_ph_diff.items()))
        
        print(hdist_to_ph_diff, "\n\n")
        ax2.plot(x, y, label=f"{bp}", marker="x", color=f"C{bp-1}")

    # ax2.set_xticks([5, 9], labels=["5", "9"])

    ax2.set_xlabel("No. sites with identical\nbase-pairing partners ", size=18)
    ax2.set_ylabel("Probability that two genotypes\nhave different phenotypes", size=18)

    for ax in (ax1, ax2):
        ax.tick_params(axis='y', which='major', labelsize=18)
        ax.tick_params(axis='y', which='minor', labelsize=18)
        ax.tick_params(axis='x', which='major', labelsize=18)
        ax.tick_params(axis='x', which='minor', labelsize=18)

    
    ax2.tick_params(axis='y', which='major', labelsize=18)
    ax2.tick_params(axis='y', which='minor', labelsize=18)
    ax2.tick_params(axis='x', which='major', labelsize=18)
    ax2.tick_params(axis='x', which='minor', labelsize=18)
    
    ax2.legend(loc="upper right", frameon=False, fancybox=False, ncol=2)

    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)
