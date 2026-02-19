#!/usr/bin/env python

import pickle
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import distance
from rna_folding.parsing import many_to_one_map_from_file_to_dict
from rna_folding.base_pairing import BasePairing
from matplotlib.ticker import MaxNLocator

def get_nc_to_gt(nc_to_gt_path):
    nc_to_gt = {}
    with open(nc_to_gt_path, "r") as f:
        for line_ in f:
            line = line_.strip().split(" ")
            nc = line[0]
            gts = line[2:]
            nc_to_gt[nc] = gts
    return nc_to_gt

def hamming_dist(gt1, gt2, l):
    d = round(distance.hamming(list(gt1), list(gt2)) * l)
    return d
    
def identical_sites(gt1, gt2, l):
    d = l-round(distance.hamming(list(gt1), list(gt2)) * l)
    return d

def get_pairwise_hdists(sample, l):
    hdist_matrix = np.zeros((len(sample), len(sample)), dtype=int)
    for i, gt1 in enumerate(sample):
        for j in np.arange(i+1, len(sample)):
            gt2 = sample[j]
            hdist_matrix[i, j] = hamming_dist(gt1, gt2, l)
    return hdist_matrix

def get_pairwise_identical_sites(sample, l):
    hdist_matrix = np.zeros((len(sample), len(sample)), dtype=int)
    for i, gt1 in enumerate(sample):
        for j in np.arange(i+1, len(sample)):
            gt2 = sample[j]
            hdist_matrix[i, j] = identical_sites(gt1, gt2, l)
    return hdist_matrix
    
def ph_diff_prob(hdist_matrix, gp_map, gt_sample):
    diff_count = {}
    diff_count_n = {}
    for hdist in np.arange(1, 13):
        diff_count[hdist] = 0
        idx1, idx2 = np.where(hdist_matrix == hdist)  # get all gt pairs with given hdist
        diff_count_n[hdist] = len(idx1)
        for i, j in zip(idx1, idx2):
            gt1 = gt_sample[i]
            gt2 = gt_sample[j]
            if gp_map.map(gt1) != gp_map.map(gt2):  # count if different ph
                diff_count[hdist] += 1
        try:
            diff_count[hdist] /= len(idx1)  # divide by number of genotypes with this hdist
        except ZeroDivisionError:
            pass
    return diff_count, diff_count_n

red_map = {1: {'J': 'J', 'K': 'R', 'L': 'R', 'M': 'M'},
               2: {'J': 'R', 'K': 'R', 'L': 'L', 'M': 'M'},
               6: {'M': 'M', 'K': 'R', 'L': 'R', 'J': 'R'},
               7: {'J': 'X', 'K': 'X', 'L': 'Y', 'M': 'Y'},  # X and why are two different redundants
               9: {'J': 'R', 'K': 'R', 'L': 'L', 'M': 'M'}}

# compute num of redundant sites and count non-redundant different sites
def map_redund_sites_to_same_letter(gts, bp_rule):
    new_gts = []
    for gt in gts:
        new_gt = ""
        for l in gt:
            l_ = red_map[bp_rule][l]
            new_gt += l_
        new_gts.append(new_gt)
    return new_gts

def get_non_neutral_neighbors(gp_map, gt):
    non_neu_neighs = []
    ph = gp_map.map(gt)
    for neigh in gp_map._neighbors(gt):
        if gp_map.map(neigh) != ph:
            non_neu_neighs.append(neigh)
    return non_neu_neighs


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

def get_avg_identical_sites_full_and_red(gp_map, nc_to_gt, nc_i,  bp, sample_size, unique_n, rng):
    nc_gts = nc_to_gt[nc_i]
    print("NC size: ", len(nc_gts))

    non_neu_neighs = get_n_non_neutral_neighbors(nc_gts, gp_map, sample_size, unique_n, rng)

    id_mtx_nc = get_pairwise_identical_sites(non_neu_neighs, l=12)
    
    if bp in red_map:
        non_neut_neighs_red = map_redund_sites_to_same_letter(non_neu_neighs, bp)
        id_mtx_nc_red = get_pairwise_identical_sites(non_neut_neighs_red, l=12)
        print(np.mean(id_mtx_nc_red.flatten()))
    else:
        id_mtx_nc_red = " "

    return non_neu_neighs, id_mtx_nc, id_mtx_nc_red


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
        gp_map = pickle.load(open(gp_map_f, "rb"))
        nc_to_gt = get_nc_to_gt(nc_to_gt_f)

        ncs = list(nc_to_gt.keys())
        rng.shuffle(ncs)

        nc_id = None
        for nc in ncs:
            if len(nc_to_gt[nc]) < 5500 and len(nc_to_gt[nc]) > 4500:
                nc_id = nc
                break
        if not nc_id:
            print("No neutral component found")
            break

        nc_gts, id_mtx, id_mtx_red = get_avg_identical_sites_full_and_red(gp_map, nc_to_gt, nc_id, bp, 1000, 200, rng)
        
        if bp in red_map:
            id_sit = id_mtx_red.flatten()[np.flatnonzero(id_mtx_red)]
            id_sit_seq = id_mtx.flatten()[np.flatnonzero(id_mtx)]
            id_sit_seq_mean = np.mean(id_sit_seq)
            def_sit_diff_nuc = [al - s for s, al in zip(id_sit_seq, id_sit)]
            # deg_sites = np.mean(id_sit)-id_sit_seq_mean
            deg_sites_mean = np.mean(def_sit_diff_nuc)
            print("ER", np.std(id_sit_seq))
            if gave_label:
                ax1.bar(bp, id_sit_seq_mean, color="teal")
                ax1.bar(bp, deg_sites_mean, bottom=id_sit_seq_mean, color="mediumaquamarine")
                
            else:
                ax1.bar(bp, deg_sites_mean, bottom=id_sit_seq_mean, color="mediumaquamarine", label="diff. nucleotides")
                ax1.bar(bp, id_sit_seq_mean, color="teal", label="ident. nucleotides")
                gave_label = True

            ax1.errorbar(bp, id_sit_seq_mean, yerr=np.std(id_sit_seq), linewidth=1, color="black", zorder=10)
            ax1.errorbar(bp + 1/10, deg_sites_mean+id_sit_seq_mean, yerr=np.std(def_sit_diff_nuc), linewidth=1, color="black", zorder=10)

            id_sites_means.append(id_sit_seq_mean)
            deg_sites_means.append(deg_sites_mean)
        else:
            id_sit = id_mtx.flatten()[np.flatnonzero(id_mtx)]
            id_sites = np.mean(id_sit)
            ax1.bar(bp, id_sites, color="teal")
            id_sites_means.append(id_sites)
            deg_sites_means.append(0)
            ax1.errorbar(bp, id_sites, yerr=np.std(id_sit), linewidth=1, color="black")
            
        diff_ph_p, diff_ph_n = ph_diff_prob(id_mtx, gp_map, nc_gts)
        y = list(diff_ph_p.values())
        x = list(diff_ph_p.keys())
        ax2.plot([x[4], x[8]], [y[4], y[8]], label=f"{bp}", marker="x", color=f"C{bp-1}")

        xs.append(" ".join([str(num) for num in x]))
        ys.append(" ".join([str(num) for num in y]))

    with open("p_over_id_sites.txt", "w") as f:
        for xstr, ystr in zip(xs, ys):
            f.write(xstr)
            f.write("\n")
            f.write(ystr)
            f.write("\n\n")

    with open("id_and_deg_sites_per_gpmap.txt", "w") as f:
        f.write(" ".join([str(num) for num in id_sites_means]))
        f.write("\n")
        f.write(" ".join([str(num) for num in deg_sites_means]))

    ax1.set_xlabel("GP map", size=18)
    ax1.set_ylabel("No. of sites with identical\nbase-pairing partners", size=18)

    ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax1.set_xticks(list(range(1, 11)))
    p =[1, 2, 3, "can.", 5, 6, 7, 8, 9, 10]
    ax1.set_xticklabels(p) #(list(range(1, 11)))

    ax2.set_xticks([5, 9], labels=["5", "9"])

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
    ax1.legend(loc="upper left", frameon=False, fancybox=False)

    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)
