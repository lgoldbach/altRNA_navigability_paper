#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pickle
import matplotlib.lines as mlines
import matplotlib as mpl
from scipy.stats.stats import pearsonr

from rna_folding.parsing import load_phenotype_and_metric_from_file, read_ruggedness_per_ph_file
from rna_folding.utils import count_bp

ph_bias = {1: .049, 2: .034, 3: .035, 4: .037, 5: .032, 6: .033, 7: .062, 8: .059, 9: .082, 10: .216}

if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--phenotype_distributions", help="Phenotype distribution files ", 
                        required=True, nargs='+')
    parser.add_argument("-o", "--output", help="Output file name "
                        "(should end in .pdf)", required=True)
    
    args = parser.parse_args()

    fig, ax = plt.subplots()

    ax.set_box_aspect(1)

    for bp, ph_dist_f in enumerate(args.phenotype_distributions, start=1):
        phenotypes, ph_count = load_phenotype_and_metric_from_file(ph_dist_f)
        ### mut group vs number of phenotypes in the top X, X and X percent.
        ph_counts_sort, ph_sort = zip(*sorted(zip(ph_count, phenotypes), reverse=True))
        ph_counts_sort = ph_counts_sort[1:]
        ph_sort = ph_sort[1:]

        folded_sum = sum(ph_counts_sort)
        top90_ph_count = 0
        top75_ph_count = 0
        top50_ph_count = 0

        s = 0
        for count in ph_counts_sort:
            s += count
            frac = s/folded_sum
            if frac < .9:
                top90_ph_count += 1
            if frac < .75:
                top75_ph_count += 1
            if frac < .5:
                top50_ph_count += 1

        # if i == 5:
        #     shift = 0.1
        # elif i == 2:
        #     shift = -0.1
        # elif i == 7:
        #     shift = -0.1
        # elif i == 4:
        #     shift = 0.1
        # else:
        #     shift = 0
        # i_ = new_order[i]
        shift = 0
        ax.scatter(ph_bias[bp], [top90_ph_count], marker="v", color=f"C{bp-1}")
        

        ax.scatter(ph_bias[bp], [top75_ph_count], color=f"C{bp-1}", marker="o")
    
        ax.scatter(ph_bias[bp], [top50_ph_count], color=f"C{bp-1}", marker="s")

        ax.vlines(x=ph_bias[bp], ymin=top90_ph_count, ymax=top50_ph_count, color=f"C{bp-1}")
        
        # if label=="natural":
        #     vl = axes[0][1].vlines(x=mut_graph_s[i_]+shift, ymin=top90_ph_count, ymax=top50_ph_count, color=col, label = "nat.")
        # else:
        #     vl = axes[0][1].vlines(x=mut_graph_s[i_]+shift, ymin=top90_ph_count, ymax=top50_ph_count, color=col, label = label)
        # ph_bias_over_mut_handles.append(vl)

    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)
