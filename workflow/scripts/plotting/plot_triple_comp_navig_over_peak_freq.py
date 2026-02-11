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

    fig, ax = plt.subplots()
    ax.set_box_aspect(1)
    
    labelsize=18
    axislabel_size=18

    colors = ["C3", "C5", "C8"]
    
    for c, navig_f, ph_dist_f in zip(colors, args.navig, args.ph_dist):
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

        ax.scatter(np.log10(x), y, color=c, alpha=.8, s=35)
        
    ax.set_xlabel("Peak phenotype frequency (log10)", fontsize=axislabel_size)
    ax.set_ylabel("Average navigability", fontsize=axislabel_size)

    ax.tick_params(axis='both', which='major', labelsize=labelsize)
    ax.tick_params(axis='both', which='minor', labelsize=labelsize)

    # ax3.legend(loc="upper left", frameon=False, fancybox=False, prop={'size': 18})
    plt.tight_layout()        
    plt.savefig(args.output, format="pdf", dpi=30)
