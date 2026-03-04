#!/usr/bin/env python

import pickle
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import distance
from rna_gpf.parsing import many_to_one_map_from_file_to_dict
from rna_gpf.base_pairing import BasePairing
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

def ph_diff_prob_from_gt_sample(gp_map, gt_sample):
    pair_count = 0
    diff_count = 0
    for i, gt1 in enumerate(gt_sample):
        for j in np.arange(i+1, len(gt_sample)):
            pair_count += 1
            gt2 = gt_sample[j]
            if gp_map.map(gt1) != gp_map.map(gt2):  # count if different ph
                diff_count += 1
    diff_frac = diff_count/pair_count
    
    # print(pair_count)
    return diff_frac


def ph_diff_prob(hdist_matrix, gp_map, gt_sample):
    diff_count = {}
    diff_count_n = {}
    for hdist in np.arange(1, 13):
        
        idx1, idx2 = np.where(hdist_matrix == hdist)  # get all gt pairs with given hdist
        diff_count_n[hdist] = len(idx1)
        if diff_count_n[hdist] < 50:
            diff_count[hdist] = np.nan  
            # print(f"not enough genotypes for distance {hdist}")
            continue
        else:
            diff_count[hdist] = 0

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


def get_n_non_neutral_neighbors(gts, gp_map, sample_size, unique_n, rng):
    non_neut_neighs = []
    for i in range(sample_size):
        nc_gt = rng.choice(gts)
        neighs = get_non_neutral_neighbors(gp_map, nc_gt)
        if neighs:
            non_neut_neigh = rng.choice(neighs)
            non_neut_neighs.append(non_neut_neigh)
    non_neut_neighs = np.unique(non_neut_neighs)[:unique_n]
    # print(f"Found: {len(non_neut_neighs)} non-neutral neighbors")
    return non_neut_neighs

def get_avg_identical_sites_full_and_red(gp_map, nc_to_gt, nc_i,  bp, sample_size, unique_n, rng):
    nc_gts = nc_to_gt[nc_i]
    # print("NC size: ", len(nc_gts))

    non_neu_neighs = get_n_non_neutral_neighbors(nc_gts, gp_map, sample_size, unique_n, rng)

    id_mtx_nc = get_pairwise_identical_sites(non_neu_neighs, l=12)
    
    if bp in red_map:
        non_neut_neighs_red = map_redund_sites_to_same_letter(non_neu_neighs, bp)
        id_mtx_nc_red = get_pairwise_identical_sites(non_neut_neighs_red, l=12)
    else:
        id_mtx_nc_red = " "

    return non_neu_neighs, id_mtx_nc, id_mtx_nc_red


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gp_map", help="GP map pickle files ", 
                        required=True, nargs='+')
    parser.add_argument("-n", "--nc_to_gt", help="neutral network to genotype map (.txt) ", 
                        required=True, nargs='+'),
    parser.add_argument("-c", "--nc_graphs", help="nc graph files", required=True, nargs="+")
    parser.add_argument("-o", "--output", help="Output file name "
                        "(should end in .pdf)", required=True)

    args = parser.parse_args()
    # fig, (ax1, ax2) = plt.subplots(ncols=2, nrows=1, figsize=(10, 5))
    fig, axes = plt.subplots(ncols=12, nrows=1, figsize=(60, 5))
    ax1 = axes[0]
    ax2 = axes[1]

    for ax in axes[1:-1]:
        # ax.grid()
        ax.set_box_aspect(1)
    
    ax1.set_box_aspect(1)
    ax2.set_box_aspect(1)

    seed = np.random.randint(low=1, high=10000)
    # seed=
    rng = np.random.default_rng(seed)
    print(f"RANDOM SEED: {seed}")

    sample_size = 1500  # number of nc genotypes for which I get non-neutral neighbors
    unique_n = 300  # final number of non-neutral neighbors picked
    nc_sample_size = 35 # number of NCs to sample

    gave_label = False

    neigh_divs_all = []
    id_sites_all = []
    for bp, (gp_map_f, nc_to_gt_f, nc_graph_f) in enumerate(zip(args.gp_map, args.nc_to_gt, args.nc_graphs), start=1):
        print(f"GP MAP {bp}")
        gp_map = pickle.load(open(gp_map_f, "rb"))
        nc_to_gt = get_nc_to_gt(nc_to_gt_f)

        ncs = list(nc_to_gt.keys())
        rng.shuffle(ncs)

        nc_ids = rng.choice(ncs, size=nc_sample_size, replace=False)

        nc_graph = pickle.load(open(nc_graph_f, "rb"))

        for nc in nc_to_gt:
            if len(nc_to_gt[nc]) != nc_graph.nodes[int(nc)]["size"]:
                print(f"NC {nc} has size {len(nc_to_gt[nc])} and {nc_graph.nodes[int(nc)]['size']}")

        
        id_sites = []  # sites with identical nucleotides
        deg_sites = []  # degenerate sites with different nucleotides
        deg_sites_all = []  # "degenerate" sites with different or identical nucleotides

        diff_ph_dict = {d: [] for d in range(1, 13)}

        avg_ph_diff_p_per_nc = []

        neigh_divs = []

        for nc_id in nc_ids:
            nc_gts, id_mtx, id_mtx_red = get_avg_identical_sites_full_and_red(gp_map, nc_to_gt, nc_id, bp, sample_size, unique_n, rng)
            
            if bp in red_map:
                id_sit_red = id_mtx_red.flatten()[np.flatnonzero(id_mtx_red)]
                id_sit_seq = id_mtx.flatten()[np.flatnonzero(id_mtx)]
                # id_sit_seq_mean = np.mean(id_sit_seq)
                id_sites += list(id_sit_seq)  # append to single superNC list
                deg_sites_all += list(id_sit_red)
                def_sit_diff_nuc = [al - s for s, al in zip(id_sit_seq, id_sit_red)]

                # deg_sites_mean = np.mean(def_sit_diff_nuc)
                deg_sites += def_sit_diff_nuc  # append to single superNC list
 
            else:
                id_sit_seq = id_mtx.flatten()[np.flatnonzero(id_mtx)]
                id_sites += list(id_sit_seq)

            # if bp in red_map:  # consider also redundant sites
            #     diff_ph_p, diff_ph_n = ph_diff_prob(id_mtx_red, gp_map, nc_gts)

            #     for d in diff_ph_p:
            #         diff_ph_dict[d].append(diff_ph_p[d])
            #     print(diff_ph_n)
            # else:
            #     diff_ph_p, diff_ph_n = ph_diff_prob(id_mtx, gp_map, nc_gts)
            #     for d in diff_ph_p:
            #         diff_ph_dict[d].append(diff_ph_p[d])
            #     print(diff_ph_n)
            
            # # compute the average phenotype diff prob, regardless of identical sites
            # avg_ph_diff_p = 0
            # n = 0
            # for d in diff_ph_n:
            #     if diff_ph_n[d] > 50:
            #         n += diff_ph_n[d]
            #         avg_ph_diff_p += (diff_ph_p[d]*diff_ph_n[d])
            # avg_ph_diff_p /= n
            # print("AAA", avg_ph_diff_p)
            # avg_ph_diff_p_per_nc.append(avg_ph_diff_p)

            diff_frac = ph_diff_prob_from_gt_sample(gp_map, nc_gts)
            # print(diff_frac)

            avg_ph_diff_p_per_nc.append(diff_frac)

            neigh_div = len(np.unique([nc_graph.nodes[neigh]["phenotype"] for neigh in nc_graph.neighbors(int(nc_id))]))
            neigh_divs.append(neigh_div)

        id_sites_mean = np.mean(id_sites)
        ax1.errorbar(bp, id_sites_mean, yerr=np.std(id_sites), linewidth=1, color="black", zorder=10)

        # print("AVG. ph. diff prob:", np.mean(avg_ph_diff_p_per_nc), avg_ph_diff_p_per_nc, id_sites_mean) 

        if bp in red_map: 
            deg_sites_mean = np.mean(deg_sites) 
            id_sites_all.append(deg_sites)
            if gave_label:
                ax1.bar(bp, id_sites_mean, color="teal")
                ax1.bar(bp, deg_sites_mean, bottom=id_sites_mean, color="mediumaquamarine")
            else:
                ax1.bar(bp, id_sites_mean, color="teal", label="ident. nucleotides")
                ax1.bar(bp, deg_sites_mean, bottom=id_sites_mean, color="mediumaquamarine", label="diff. nucleotides")
                gave_label = True
            ax1.errorbar(bp + 1/10, deg_sites_mean+id_sites_mean, yerr=np.std(deg_sites), linewidth=1, color="black", zorder=10)

            axes[-1].errorbar(np.mean(deg_sites_all), np.mean(avg_ph_diff_p_per_nc), yerr=np.std(avg_ph_diff_p_per_nc), xerr=np.std(deg_sites_all), label=f"{bp}", color=f"C{bp-1}", fmt='o')
            print(" x ", np.mean(deg_sites_all), " y ", np.mean(avg_ph_diff_p_per_nc), " yerr ", np.std(avg_ph_diff_p_per_nc), " x ",np.std(deg_sites_all))

            axes[-2].errorbar(np.mean(deg_sites_all), np.mean(neigh_divs), yerr=np.std(neigh_divs), xerr=np.std(deg_sites_all), label=f"{bp}", color=f"C{bp-1}", fmt='o')
        else:
            id_sites_all.append(id_sites)
            if gave_label:
                ax1.bar(bp, id_sites_mean, color="teal")
            else:
                ax1.bar(bp, id_sites_mean, color="teal", label="ident. nucleotides")
                gave_label = True
            axes[-1].errorbar(id_sites_mean, np.mean(avg_ph_diff_p_per_nc), yerr=np.std(avg_ph_diff_p_per_nc), xerr=np.std(id_sites), label=f"{bp}", color=f"C{bp-1}", fmt='o')
            axes[-2].errorbar(id_sites_mean, np.mean(neigh_divs), yerr=np.std(neigh_divs), xerr=np.std(id_sites), label=f"{bp}", color=f"C{bp-1}", fmt='o')
            
            print(" x ", id_sites_mean, "y ", np.mean(avg_ph_diff_p_per_nc), " yerr ", np.std(avg_ph_diff_p_per_nc), " xerr ", np.std(id_sites))
        print("neigh divs mean and std", np.mean(neigh_divs), np.std(neigh_divs))
        neigh_divs_all.append(neigh_divs)
        
        
    
    with open("neigh_divs_from_NC_sample_fig6.txt", "w") as f:
        for divs in neigh_divs_all:
            f.write(" ".join([str(d) for d in divs]))
            f.write("\n")

    with open("id_sites_from_NC_sample_fig6.txt", "w") as f:
        for id_ in id_sites_all:
            f.write(" ".join([str(d) for d in id_sites_all]))
            f.write("\n")

        # x = []
        # y = []
        # err = []
        # for d in diff_ph_dict:
        #     x.append(d)
        #     diff_p_wo_nan = [v for v in diff_ph_dict[d] if v != np.nan]
        #     print(d, diff_p_wo_nan, np.mean(diff_p_wo_nan), np.nanmean(diff_p_wo_nan), np.nanmean(diff_ph_dict[d]))
        #     y.append(np.nanmean(diff_ph_dict[d]))

        #     err.append(np.nanstd(diff_ph_dict[d]))
        #     # if diff_p_wo_nan:
        #     #     print("Passed")
        #     #     y.append(np.mean(diff_p_wo_nan))
        #     #     err.append(np.std(diff_p_wo_nan))
        #     # else:
        #     #     print("Not passed")
        #     #     y.append(0)
        #     #     err.append(0)

        # axes[bp].errorbar(x, y, yerr=err, label=f"{bp}", marker="x", color=f"C{bp-1}")
        # # ax2.plot([x[4], x[8]], [y[4], y[8]], label=f"{bp}", marker="x", color=f"C{bp-1}")
        # ax2.plot(x, y, label=f"{bp}", marker="x", color=f"C{bp-1}")


    ax1.set_xlabel("GP map", size=18)
    ax1.set_ylabel("No. of sites with identical\nbase-pairing partners", size=18)

    ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax1.set_xticks(list(range(1, 11)))
    p =[1, 2, 3, "can.", 5, 6, 7, 8, 9, 10]
    ax1.set_xticklabels(p) #(list(range(1, 11)))

    # ax2.set_xticks([5, 9], labels=["5", "9"])

    for ax in (ax1, ax2, axes[-1], axes[-2]):
        ax.tick_params(axis='y', which='major', labelsize=18)
        ax.tick_params(axis='y', which='minor', labelsize=18)
        ax.tick_params(axis='x', which='major', labelsize=18)
        ax.tick_params(axis='x', which='minor', labelsize=18)

    ax1.legend(loc="upper left", frameon=False, fancybox=False)

    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)
