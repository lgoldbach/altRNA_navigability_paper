#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pickle
import matplotlib.lines as mlines
import matplotlib as mpl
from scipy.stats.stats import pearsonr

from rna_gpf.parsing import load_phenotype_and_metric_from_file, read_ruggedness_per_ph_file
from rna_gpf.utils import count_bp
from rna_gpf.analysis import get_peaks


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--peak_sizes", help="Peak sizes files ", 
                        required=True, nargs='+')
    parser.add_argument("-n", "--nc_graphs", help="nc graph pickle files ", 
                        required=True, nargs='+')
    parser.add_argument("-o", "--output", help="Output file name "
                        "(should end in .pdf)", required=True)

    args = parser.parse_args()
    fig, ax = plt.subplots()
    
    ax.set_box_aspect(1)

    rng = np.random.default_rng()
    for i, (peak_sizes_f, nc_graph_f) in enumerate(zip(args.peak_sizes, args.nc_graphs)):
        glo_peak_ph_per_fl, peak_sizes_per_fl = load_phenotype_and_metric_from_file(peak_sizes_f, ignore="............", dtype=int)
        nc_graph = pickle.load(open(nc_graph_f, "rb"))

        sum_of_folded = sum([nc_graph.nodes[node]["size"] for node in nc_graph])

        ## below is code for generating 1000 random fitness landscapes and computing average
        phenotypes_peak = np.unique(glo_peak_ph_per_fl)
        phenotypes = np.unique([nc_graph.nodes[node]["phenotype"] for node in nc_graph])
        all_peak_sizes = []
        low_f = 0
        high_f = 1
        for i in range(10000):
            ph_to_f = {}
            # randomly assign fitness from the open interval (low_f, upp_f), i.e. for
            for ph in phenotypes:
                f = np.random.uniform(low_f, high_f)
                ph_to_f[ph] = f

            peaks_nc, peaks_f = get_peaks(nc_graph, ph_to_f)
            peak_size = 0
            for nc in peaks_nc:
                peak_size += nc_graph.nodes[nc]["size"]
            all_peak_sizes.append(peak_size)
        peak_size_avg_numeric = np.mean(all_peak_sizes)  / sum_of_folded

        ### compute analytical prediction
        peak_size_pred = 0
        for node in nc_graph.nodes:
            nc_size = nc_graph.nodes[node]["size"]
            evo = len(np.unique([nc_graph.nodes[neigh]["phenotype"] for neigh in nc_graph.neighbors(node)]))
            peak_size_pred += nc_size/(evo+1)

        peak_size_pred /= sum_of_folded

        # peak_size_avg_numeric = np.mean([peak_size/sum_of_folded for peak_size in peak_sizes_per_fl])

        ax.scatter(peak_size_pred, peak_size_avg_numeric)
        ax.set_xlabel("Global + local peak fraction\n(analytical prediction)", size=15)
        ax.set_ylabel("Global + local peak fraction\n(numerically estimate)", size=15)
    
        ax.plot([0, 1], [0, 1], transform=ax.transAxes, color="0.2", zorder=-10, linewidth=0.5, linestyle="--")
        ax.tick_params(axis='y', which='major', labelsize=15)
        ax.tick_params(axis='y', which='minor', labelsize=15)
        ax.tick_params(axis='x', which='major', labelsize=15)
        ax.tick_params(axis='x', which='minor', labelsize=15)
        
    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)
