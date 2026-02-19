#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

from rna_gpf.parsing import load_phenotype_and_metric_from_file, read_navigability_per_fl
from rna_gpf.utils import ranked_ph_distribution


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--phenotype_dist", help="Path to phenotype "
                        "distribution files", required=True)
    parser.add_argument("-r", "--ref", help="Reference phenotype distribution  "
                        "files ", required=True)
    parser.add_argument("-o", "--output", help="Output file name "
                        "(should end in .pdf)", required=True)
    parser.add_argument("-n", "--nc", help="Path to neutral component file", required=True)
    parser.add_argument("-m", "--ref_nc", help="Reference neutral component "
                        "file ", required=True)
    # parser.add_argument("-c", "--rank_cutoff", help="How many ranks to plot max.",
                        # type=int, required=False)
    parser.add_argument("--walks", help="Walk success of query", 
                        required=True, 
                        type=str)
    parser.add_argument("--refwalks", help="Walk success of query", 
                        required=True, 
                        type=str)
    parser.add_argument("--sample_size", help="How big is one set of samples per phenotype, i.e. how many samples are takend for a given fitness landscape instance of a phenotype",
                        required=True, type=int)


    args = parser.parse_args()

    fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(15, 5))
    ax1.set_box_aspect(1)
    ax2.set_box_aspect(1)
    ax3.set_box_aspect(1)
    labelsize=18
    axislabel_size=18
    ### phenotype bias part
    phenotypes, distr = ranked_ph_distribution(ph_distr_file=args.ref,
                                                log=True)

    distr = distr[1:]  # remove unfolded phenotype frequency
    x = range(1, distr.shape[0]+1)
    
    # sc = ax1.plot(x, distr, marker="", label=f"ViennaRNA", color="black", linewidth=7)
    sc = ax1.scatter(x, distr, color="0.2", zorder=10,  label=f"ViennaRNA g-p map", marker="s", s=30) # , edgecolor="black", linewidth=0.5)

    phenotypes, distr = ranked_ph_distribution(ph_distr_file=args.phenotype_dist,
                                                log=True)
    distr = distr[1:]  # remove unfolded phenotype frequency
    x = range(1, distr.shape[0]+ 1)
    
    # sc = ax1.plot(x, distr, marker="", label=f"Global ranking", color="0.4", linewidth=7, linestyle="dotted", zorder=10)
    sc = ax1.scatter(x, distr, color="orange", zorder=10, label=f"mfe-score g-p map", s=35) # edgecolor="black", linewidth=0.5)

    ax1.tick_params(axis='both', which='major', labelsize=labelsize)
    ax1.tick_params(axis='both', which='minor', labelsize=labelsize)

    ax1.set_xlabel("Rank", fontsize=axislabel_size)
    ax1.set_ylabel("Phenotype frequency (log10)", fontsize=axislabel_size)
    # ax1.grid()
    # ax1.legend(loc="lower left", prop={'size': 18}, frameon=False)

    ### neutral component part
    nc_sizes = []
    with open(args.ref_nc, "r") as file:
        for line in file:
            for size in line.strip().split(" "):
                if size:
                    nc_sizes.append(int(size))

    nc_sizes_sort = np.sort(nc_sizes)[::-1]
    nc_sizes = nc_sizes_sort

    # if args.rank_cutoff:
    #     try:
    #         nc_sizes = nc_sizes_sort[:args.rank_cutoff]
    #     except IndexError:
    #         pass
    
    x = range(nc_sizes.shape[0])
    # ax2.plot(x[1:], np.log10(nc_sizes)[1:], marker="", label=f"ViennaRNA", color="black", linewidth=7, markersize=5)
    sc = ax2.scatter(x[1:], np.log10(nc_sizes)[1:], color="0.2", zorder=10,  label=f"ViennaRNA g-p map", marker="s", s=30)

    nc_sizes = []
    with open(args.nc, "r") as file:
        nc_sizes = []
        for line in file:
            for size in line.strip().split(" "):
                if size:
                    nc_sizes.append(int(size))

    nc_sizes_sort = np.sort(nc_sizes)[::-1]
    nc_sizes = nc_sizes_sort
    # if args.rank_cutoff:
    #     try:
    #         nc_sizes = nc_sizes_sort[:args.rank_cutoff]
    #     except IndexError:
    #         pass

    step = 20
    nc_sizes = nc_sizes[1:]

    step_nc_sizes = [nc_sizes[i] for i in range(0, len(nc_sizes), step)]  
    
    step_x = range(0, len(nc_sizes), step)
    # ax2.plot(step_x, np.log10(step_nc_sizes), marker="", label=f"Global ranking", color="0.4", linewidth=7, zorder=10, linestyle="dotted")
    sc = ax2.scatter(np.arange(start=1, stop=len(nc_sizes)+1), np.log10(nc_sizes), color="orange", zorder=10, label=f"mfe-score g-p map", s=35)
    ax2.tick_params(axis='both', which='major', labelsize=labelsize)
    ax2.tick_params(axis='both', which='minor', labelsize=labelsize)

    ax2.set_xlabel("Rank", fontsize=axislabel_size)
    ax2.set_ylabel("Neutral component size (log10)", fontsize=axislabel_size)
    # ax2.grid()
    # ax2.legend(loc="lower left", prop={'size': 18}, frameon=False)

    
    ### navigability over frequency part
    ref_ph, ref_freq = load_phenotype_and_metric_from_file(args.ref, dtype=int, ignore="............")
    ph, freq = load_phenotype_and_metric_from_file(args.phenotype_dist, dtype=int, ignore="............")
    
    sum_f = np.sum(ref_freq)
    ref_freq = [fre/sum_f for fre in ref_freq]
    sum_f = np.sum(freq)
    freq = [fre/sum_f for fre in freq]

    ref_freq_sort, ref_ph_sort = zip(*sorted(zip(ref_freq, ref_ph)))  # sort both lists by frequency
    freq_sort, ph_sort = zip(*sorted(zip(freq, ph)))  # sort both lists by frequency

    ph_sort = np.array(ph_sort)
    ref_ph_sort = np.array(ref_ph_sort)

    navigs_per_ph = read_navigability_per_fl(args.walks)
    ref_navigs_per_ph = read_navigability_per_fl(args.refwalks)


    ref_l = "ViennaRNA"
    query_l = "Sequence-free method"

    ref_x = []
    ref_y = []

    for ph, f in zip(ref_ph_sort, ref_freq_sort):
        if ph in ref_navigs_per_ph:
            ref_y.append(np.mean(ref_navigs_per_ph[ph]))
            ref_x.append(f)
    
    ax3.scatter(np.log10(ref_x), ref_y, label=ref_l, marker='s', color="0.2", s=30)

    x = []
    y = []

    for ph, f in zip(ph_sort, freq_sort):
        if ph in navigs_per_ph:
            y.append(np.mean(navigs_per_ph[ph]))
            x.append(f)

    
    ax3.scatter(np.log10(x), y, label=query_l, color="orange", alpha=.8, s=35)
        
    
    # x_query = []
    # y_query = []
    # y_err = []
    # y_err_2d = [[],[]]
    # for p in ph_sort:
    #     if p in walk_success and p in query_d:
    #         s = walk_success[p]
    #         m = np.median(s)
    #         p_low = np.abs(m-np.percentile(s, q=25))
    #         p_high = np.abs(m-np.percentile(s, q=75))

    #         y_query.append(m)
    #         x_query.append(query_d[p])
    #         y_err.append(np.std(s))

    #         y_err_2d[0].append(p_low)
    #         y_err_2d[1].append(p_high)
            
    
    # ref_l = "ViennaRNA"
    # query_l = "Global ranking"
    
    # # ax3.errorbar(np.log10(x_ref), np.array(y_ref), yerr=np.array(y_err_2d_ref), label=ref_l, linestyle='', marker='s', elinewidth=0.2, color="0.2", alpha=1, markeredgewidth=0, markersize=7)
    # # ax3.errorbar(np.log10(x_query), np.array(y_query), yerr=np.array(y_err_2d), label=query_l, linestyle='', marker='o', elinewidth=0.2, color="orange", alpha=0.8, markeredgewidth=0, markersize=7)
    # print(y_ref, y_query)
    # ax3.errorbar(np.log10(x_ref), np.array(y_ref), yerr=0, label=ref_l, linestyle='', marker='s', elinewidth=0.2, color="0.2", alpha=1, markeredgewidth=0, markersize=7)
    # ax3.errorbar(np.log10(x_query), np.array(y_query), yerr=0, label=query_l, linestyle='', marker='o', elinewidth=0.2, color="orange", alpha=0.8, markeredgewidth=0, markersize=7)
  

    ax3.set_xlabel("Peak phenotype frequency (log10)", fontsize=axislabel_size)
    ax3.set_ylabel("Average navigability", fontsize=axislabel_size)


    ax3.tick_params(axis='both', which='major', labelsize=labelsize)
    ax3.tick_params(axis='both', which='minor', labelsize=labelsize)

    # ax3.legend(loc="upper left", frameon=False, fancybox=False, prop={'size': 18})
    plt.tight_layout()        
    plt.savefig(args.output, format="pdf", dpi=30)
