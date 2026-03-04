#!/usr/bin/env python

import pickle
import argparse
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.ticker import MaxNLocator

degen = {1: .17, 2: .17, 3: 0, 4: 0, 5: 0, 6: .5, 7: .33, 8: 0, 9: .17, 
            10: 0}

if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--nc_graphs", help="NC graph pickle files ", 
                        required=True, nargs='+')
    parser.add_argument("-o", "--output", help="Output file name "
                        "(should end in .pdf)", required=True)

    args = parser.parse_args()
   
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(10, 5))
    
    for ax in (ax1, ax2):
        ax.set_box_aspect(1)
        ax.tick_params(axis='y', which='major', labelsize=18)
        ax.tick_params(axis='y', which='minor', labelsize=18)
        ax.tick_params(axis='x', which='major', labelsize=18)
        ax.tick_params(axis='x', which='minor', labelsize=18)

    for bp, nc_graph_f in enumerate(args.nc_graphs, start=1):
        nc_graph = pickle.load(open(nc_graph_f, "rb"))
        neigh_divs = []
        for nc in nc_graph:
            neigh_div = len(np.unique([nc_graph.nodes[neigh]["phenotype"] for neigh in nc_graph.neighbors(nc)]))
            neigh_divs.append(neigh_div)
        
        neigh_divs_mean = np.mean(neigh_divs)
        neigh_divs_std = np.std(neigh_divs)
        ax1.bar(bp, neigh_divs_mean, yerr=neigh_divs_std, color="grey")

        ax2.scatter(degen[bp], neigh_divs_mean, color=f"C{bp-1}", label=bp)
        ax2.errorbar(degen[bp], neigh_divs_mean, yerr=neigh_divs_std, color=f"C{bp-1}")

    ax1.set_xlabel("GP map", size=18)
    ax1.set_ylabel("Neighborhood diversity", size=18)

    ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax1.set_xticks(list(range(1, 11)))
    p =[1, 2, 3, "can.", 5, 6, 7, 8, 9, 10]
    ax1.set_xticklabels(p) #(list(range(1, 11)))

    ax2.legend()
    
    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)
