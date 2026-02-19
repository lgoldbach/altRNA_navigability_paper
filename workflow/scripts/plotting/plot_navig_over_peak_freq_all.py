#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

from rna_gpf.parsing import load_phenotype_and_metric_from_file, read_navigability_per_fl
from rna_gpf.utils import ranked_ph_distribution


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--ph_dist", help="Path to phenotype "
                        "distribution files", nargs="+", required=True)
    parser.add_argument("-n", "--navig", help="Files with navigability per fl", 
                        required=True, nargs="+", type=str)
    parser.add_argument("-i", "--ignore", help="phenotype to ignore", required=False, type=str, default=None)
    parser.add_argument("-o", "--output", help="Output file name "
                        "(should end in .pdf)", required=True)

    args = parser.parse_args()

    fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(15, 15))
    for ax in axes.flat:
        ax.set_box_aspect(1)

    labelsize=18
    axislabel_size=18
    
    id_mod = 0  # keep track of whether we have plotted canonical data yet to plot inside right plot
    axes_flat = axes.flat
    for i, (navig_f, ph_dist_f) in enumerate(zip(args.navig, args.ph_dist)):
        ph, freq = load_phenotype_and_metric_from_file(ph_dist_f, dtype=int, ignore=args.ignore)

        sum_f = np.sum(freq)
        freq = [fre/sum_f for fre in freq]  # normalize frequencies

        navigs_per_ph = read_navigability_per_fl(navig_f)

        x = []
        y = []

        for ph, f in zip(ph, freq):
            if ph in navigs_per_ph:
                y.append(np.mean(navigs_per_ph[ph]))
                x.append(f)
        
        c = f"C{i}"
        
        if i == 3:   # plot canonical data in all plots as reference
            for ax_ in axes_flat:
                ax_.scatter(np.log10(x), y, color=c, alpha=.8, s=35, label="canon.")
                ax_.vlines(np.log10(np.mean(x)), -0.01, 1.01, color=c, zorder=-10)
            id_mod = -1
            continue

        
        ax = axes_flat[i+id_mod]

        label = f"{i+1}"

        ax.scatter(np.log10(x), y, color=c, alpha=.8, s=35, label=label)
        ax.vlines(np.log10(np.mean(x)), -0.01, 1.01, color=c, zorder=-10)
    # ax.scatter(np.log10([10**-1.1, 10**-1.2, 10**-1.5, 10**-1.7]), [0.9]*4)
    # ax.set_xlabel("Peak phenotype frequency (log10)", fontsize=axislabel_size)
    # ax.set_ylabel("Average navigability", fontsize=axislabel_size)

        ax.set_ylim(-0.01, 1.01)
        # ax.set_xlim(-4.3, -0.2)

        

        ax.tick_params(axis='both', which='major', labelsize=labelsize)
        ax.tick_params(axis='both', which='minor', labelsize=labelsize)

    for ax in axes_flat:
        ax.legend(title="GP map", frameon=False, fancybox=False, fontsize=15, title_fontsize=17)
    # ax3.legend(loc="upper left", frameon=False, fancybox=False, prop={'size': 18})
    plt.tight_layout()        
    plt.savefig(args.output, format="pdf", dpi=30)
