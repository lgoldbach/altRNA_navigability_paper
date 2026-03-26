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
    fig, axes = plt.subplots(nrows=4, ncols=3, figsize=(15, 20))
    
    for ax in axes.flat:
        ax.set_box_aspect(1)
    
    ax = axes[0][0]
    ax_id_mod = 0

    axes_flat = axes[1:].flatten()
    rng = np.random.default_rng()

    loc_peak_size_preds = []
    loc_peak_size_avg_across_fl_vs = []
    for i, (peak_sizes_f, nc_graph_f) in enumerate(zip(args.peak_sizes, args.nc_graphs)):
        glo_peak_ph_per_fl, peak_sizes_per_fl = load_phenotype_and_metric_from_file(peak_sizes_f, ignore="............", dtype=int)
        nc_graph = pickle.load(open(nc_graph_f, "rb"))


        sum_of_folded = sum([nc_graph.nodes[node]["size"] for node in nc_graph])

        ph_to_count = {}
        nc_sizes = []
        for nc in nc_graph:
            ph = nc_graph.nodes[nc]["phenotype"]
            if ph not in ph_to_count:
                ph_to_count[ph] = nc_graph.nodes[nc]["size"]
            else:
                ph_to_count[ph] += nc_graph.nodes[nc]["size"]
            nc_sizes.append(nc_graph.nodes[nc]["size"])
        print(i+1, np.median(nc_sizes), np.median(nc_sizes)*len(nc_sizes), np.mean(nc_sizes), np.mean(nc_sizes)/len(nc_sizes), (np.mean(nc_sizes)/len(nc_sizes))/sum_of_folded)
        exp_glo_peak = np.mean(list(ph_to_count.values()))/sum_of_folded
        # print(i, exp_glo_peak)
        ## below is code for generating 1000 random fitness landscapes and computing average
        phenotypes_peak = np.unique(glo_peak_ph_per_fl)
        phenotypes = np.unique([nc_graph.nodes[node]["phenotype"] for node in nc_graph])
        print(len(phenotypes))
        all_peak_sizes = []
        loc_peak_sizes = []
        low_f = 0
        high_f = 1
        max_f = 0
        for it in range(200):
            ph_to_f = {}
            # randomly assign fitness from the open interval (low_f, upp_f), i.e. for
            for ph in phenotypes:
                f = np.random.uniform(low_f, high_f)
                ph_to_f[ph] = f
                if f > max_f:
                    max_f = f

            peaks_nc, peaks_f = get_peaks(nc_graph, ph_to_f)
            peak_size = 0
            glo_peak_size = 0
            for nc in peaks_nc:
                nc_s = nc_graph.nodes[nc]["size"]
                peak_size += nc_s
                if ph_to_f[nc_graph.nodes[nc]["phenotype"]] == max_f:
                    glo_peak_size += nc_s
            all_peak_sizes.append(peak_size)
            loc_peak_sizes.append(peak_size - glo_peak_size)

        peak_size_avg_numeric = np.mean(all_peak_sizes)  / sum_of_folded
        loc_peak_size_avg_numeric = np.mean(loc_peak_sizes)  / sum_of_folded

        ### compute analytical prediction
        peak_size_pred = 0
        loc_peak_size_pred = 0
        glo_peak_size_pred = 0
        evos = []
        nc_sizes = []
        for node in nc_graph.nodes:
            nc_size = nc_graph.nodes[node]["size"]
            evo = len(np.unique([nc_graph.nodes[neigh]["phenotype"] for neigh in nc_graph.neighbors(node)]))
            peak_size_pred += nc_size/(evo+1)   

            peak_prob = 1/(evo+1)
            glo_prob = 1/(len(phenotypes)+1)
            loc_peak_prob = peak_prob - glo_prob
            # loc_peak_size_pred += nc_size * loc_peak_prob

            loc_peak_size_pred += ((nc_size/(evo+1)) - (nc_size/len(phenotypes)))

    
            glo_peak_size_pred += nc_size * glo_prob

            evos.append(evo)
            nc_sizes.append(nc_size)
        
        nc_sizes_sort, evos_sort = zip(*sorted(zip(nc_sizes, evos)))
        peak_contrib_per_nc_sort = [(nc_size_i/(evo_i+1)) for nc_size_i, evo_i in zip(nc_sizes_sort, evos_sort)]
        cum_peak_pred = [(nc_sizes_sort[0]/(evos_sort[0]+1))/sum_of_folded]
        for nc_size, evo in zip(nc_sizes_sort[1:], evos_sort[1:]):
            cum_peak_pred.append(cum_peak_pred[-1]+((nc_size/(evo+1))/sum_of_folded))
        
        # if i == 3:
        #     ax_id_mod = -1
        #     for axi in axes_flat:  # plot canonical everywhere
        #         # axes_flat[i].scatter(np.log10(nc_sizes), [(nc_size/(evo+1))/sum_of_folded for nc_size, evo in zip(nc_sizes, evos)])
        #         # axi.plot(np.log10(nc_sizes_sort), cum_peak_pred, marker=".", color=f"C{i}")
        #         axi.plot(np.log10(nc_sizes_sort), cum_peak_pred, marker=".", color=f"C{i}", zorder=-10)
        #         # axi2 = axi.twinx()
        #         # axi2.scatter(np.log10(nc_sizes), evos)
        # else:
        #     # axes_flat[i].scatter(np.log10(nc_sizes), evos)
        #     axi = axes_flat[i+ax_id_mod]
        #     # axi.plot(np.log10(nc_sizes_sort), cum_peak_pred, marker=".", color=f"C{i}")
        #     axi.plot(np.log10(nc_sizes_sort), cum_peak_pred, marker=".", color=f"C{i}")
        #     axi.set_xlim(2, 6)
        #     axi.set_ylim(0, 0.06)
            
        #     # axi2 = axi.twinx()
        #     # axi2.scatter(np.log10(nc_sizes), evos)

        if i == 3:
            # pass
            ax_id_mod = -1
            for axi in axes_flat:  # plot canonical everywhere
                # axes_flat[i].scatter(np.log10(nc_sizes), [(nc_size/(evo+1))/sum_of_folded for nc_size, evo in zip(nc_sizes, evos)])
                # axi.plot(np.log10(nc_sizes_sort), cum_peak_pred, marker=".", color=f"C{i}")
                axi.plot(np.log10(nc_sizes_sort), np.log10(peak_contrib_per_nc_sort), marker=".", color=f"C{i}", zorder=-10)
                # axi2 = axi.twinx()
                # axi2.scatter(np.log10(nc_sizes), evos)
        else:
            # axes_flat[i].scatter(np.log10(nc_sizes), evos)
            axi = axes_flat[i+ax_id_mod]
            # axi.plot(np.log10(nc_sizes_sort), cum_peak_pred, marker=".", color=f"C{i}")
            axi.plot(np.log10(nc_sizes_sort), np.log10(peak_contrib_per_nc_sort), marker=".", color=f"C{i}")
            # axi.set_ylim(0, .06)
            # if i != 2:
            axi.set_ylim(0, 4.7)
            axi.set_xlim(2, 6)
            
            axi2 = axi.twinx()
            axi2.scatter(np.log10(nc_sizes_sort), cum_peak_pred, color=f"C{i}")
            axi2.set_ylim(0, 0.07)

        loc_peak_size_pred_by_div = peak_size_pred - (sum_of_folded/len(phenotypes))
        loc_peak_size_pred_by_div /= sum_of_folded
        peak_size_pred /= sum_of_folded
        # loc_peak_size_pred -= (sum_of_folded/len(phenotypes))
        loc_peak_size_pred /= sum_of_folded
        glo_peak_size_pred /= sum_of_folded
        

        # take average peak size from sample
        peak_size_avg_across_fl = np.mean([peak_size/sum_of_folded for peak_size in peak_sizes_per_fl])
        loc_peak_size_avg_across_fl = []
        glo_peak_size_avg_across_fl = []
        for ph, peak_size in zip(glo_peak_ph_per_fl, peak_sizes_per_fl):
            glo_peaks = sum([nc_graph.nodes[nc]["size"] for nc in nc_graph if nc_graph.nodes[nc]["phenotype"]==ph])
            glo_peak_size_avg_across_fl.append(glo_peaks)
            loc_peak_size = peak_size-glo_peaks
            loc_peak_size_avg_across_fl.append(loc_peak_size)

        glo_peak_size_avg_across_fl_v = np.mean(glo_peak_size_avg_across_fl) / sum_of_folded
        loc_peak_size_avg_across_fl_v = np.mean(loc_peak_size_avg_across_fl) / sum_of_folded
            
        loc_peak_size_preds.append(loc_peak_size_pred)
        loc_peak_size_avg_across_fl_vs.append(loc_peak_size_avg_across_fl_v)
        
        ax.scatter(loc_peak_size_pred, loc_peak_size_avg_across_fl_v)
        ax.set_xlabel("Global + local peak fraction\n(analytical prediction)", size=15)
        ax.set_ylabel("Global + local peak fraction\n(numerically estimate)", size=15)
    
        # ax.plot([0, 1], [0, 1], transform=ax.transAxes, color="0.2", zorder=-10, linewidth=0.5, linestyle="--")
        ax.plot([0.001, 0.05], [0.001, 0.05], color="0.2", zorder=-10, linewidth=0.5, linestyle="--")
        ax.tick_params(axis='y', which='major', labelsize=15)
        ax.tick_params(axis='y', which='minor', labelsize=15)
        ax.tick_params(axis='x', which='major', labelsize=15)
        ax.tick_params(axis='x', which='minor', labelsize=15)

    r, p = pearsonr(loc_peak_size_preds, loc_peak_size_avg_across_fl_vs) 
    print(r, p)
    p_str = "%.3g" % p
    ax.text(0, .9, f'r = {np.round(r, 2)}\np = {p_str}', transform=ax.transAxes)

    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)
