#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pickle
import matplotlib.lines as mlines
import matplotlib as mpl

from rna_folding.parsing import load_phenotype_and_metric_from_file, read_ruggedness_per_ph_file
from rna_folding.utils import count_bp

ph_bias = {1: .049, 2: .034, 3: .035, 4: .037, 5: .032, 6: .033, 7: .062, 8: .059, 9: .082, 10: .216}

if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--phenotype_distributions", help="Phenotype distribution files ", 
                        required=True, nargs='+')
    parser.add_argument("-n", "--nc_graphs", help="NC graphs", required=True, nargs="+")
    parser.add_argument("-o", "--output", help="Output file name "
                        "(should end in .pdf)", required=True)
    
    args = parser.parse_args()

    fig, axes = plt.subplots(nrows=5, ncols=2, figsize=(13, 27))
    axes_flat = axes.flat

    axes2 = []
    # for ax in axes_flat:
    #     ax2 = ax.twinx()
    #     axes2.append(ax2)

    labelsize = 20

    # print(len(axes_flat), len(axes2))

    # print(list(zip(axes_flat, args.phenotype_distributions, args.nc_graphs)))

    for bp, (ax, ph_dist_f, nc_graph_f) in enumerate(zip(axes_flat, args.phenotype_distributions, args.nc_graphs), start=1):
        phenotypes, ph_count = load_phenotype_and_metric_from_file(ph_dist_f)
        
        ph_counts_sort, ph_sort = zip(*sorted(zip(ph_count, phenotypes), reverse=True))
        ph_counts_sort = ph_counts_sort[1:]
        ph_sort = ph_sort[1:]

        folded_sum = sum(ph_counts_sort)
        all_gt_count = 4**12

        ph_frac_of_sort = [c/all_gt_count for c in ph_counts_sort]

        nc_graph = pickle.load(open(nc_graph_f, "rb"))
        nc_count = len(nc_graph.nodes)

        # count in how many neutral networks each phenotype is a neighbor of
        num_of_nc_ph_appears_in = {ph: 0 for ph in ph_sort}
   
        for nc in nc_graph:  # go through each NC
            neigh_phs = np.unique([nc_graph.nodes[neigh]["phenotype"] for neigh in nc_graph.neighbors(nc)])  # get uniq. list of its neighbor's ph.
            for neigh_ph in neigh_phs: # count towards that phenotypes number of NCs
                num_of_nc_ph_appears_in[neigh_ph] += 1  

        # divide by number of neutral components
        frac_of_nc_ph_appears_in = [num_of_nc_ph_appears_in[ph]/nc_count for ph in ph_sort]

        frac_of_nc_ph_appears_in_cum = []
        frac_cum = 0
        for frac in frac_of_nc_ph_appears_in:
            frac_cum += frac
            frac_of_nc_ph_appears_in_cum.append(frac_cum)

    
        ax.bar(range(len(ph_sort)), frac_of_nc_ph_appears_in, color="0.4")

        ax2 = ax.twinx()
        axes2.append(ax2)
        ax2.scatter(range(len(ph_sort)), np.log10(ph_frac_of_sort), color="tab:blue", edgecolor="white", s=40, linewidth=.5)
        # ax2.plot(range(len(ph_sort)), frac_of_nc_ph_appears_in_cum, marker="o", color="black")
        
        ax.set_box_aspect(1)
        ax.set_xlim(-1, 37)
        ax.set_ylim(0, 1.09)

        ax2.set_ylim(-6.5, -1)
        
        for ax_ in (ax, ax2):
            ax_.tick_params(axis='both', which='major', labelsize=labelsize)
            ax_.tick_params(axis='both', which='minor', labelsize=labelsize)
            ax_.set_box_aspect(1)    
        
        label = f"GP map {bp}"
        if bp == 4:
            label = "canon. GP map"

        neigh_div_avg = np.round(sum(frac_of_nc_ph_appears_in), 2)
        ax.text(x=0.98, y=0.85, ha="right", s=f"{label}\n{neigh_div_avg}", size=22, transform=ax.transAxes)

    for ax in axes[-1]:  # last row
        ax.set_xlabel("Phenotype rank\n(desc. order of frequency)", size=labelsize)

    for ax in axes[:, 0]:  # left col
        ax.set_ylabel("Frac. of neut. net. neighborhoods\n that contain phenotype", size=labelsize)

    for ax2 in axes2[1::2]:  # left col
        ax2.set_ylabel("Phenotype frequency (log10)", size=labelsize, color="tab:blue")

    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)
