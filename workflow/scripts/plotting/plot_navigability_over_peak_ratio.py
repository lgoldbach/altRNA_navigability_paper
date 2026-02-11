#!/usr/bin/env python

import argparse
import pickle
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats.stats import pearsonr

from rna_gpf.parsing import load_phenotype_and_metric_from_file, read_navigability_per_fl



if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--navigability", help="Navigability files", nargs="+", required=True)
    parser.add_argument("-p", "--peak_sizes", help="Peak size files", nargs="+", required=True)
    parser.add_argument("-f", "--ph_dist", help="phenotype distribution files", nargs="+", required=True)
    parser.add_argument("-o", "--output", help="Output file for plot (.pdf)",
                        required=True)
    
    args = parser.parse_args()
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(10, 5), sharey=False)
    ax1.set_box_aspect(1)
    ax2.set_box_aspect(1)

    mean_peak_ratios = []
    mean_navigs = []
    mean_loc_peak_fracs = []
    for i, (n, p, f) in enumerate(zip(args.navigability, args.peak_sizes, args.ph_dist), start=1):
        ph_f, ph_counts = load_phenotype_and_metric_from_file(f, dtype=int, ignore="............")
        ph_to_counts = dict(zip(ph_f, ph_counts))
        folded_count = sum(ph_counts)

        ph_n, navigs_per_fl =  load_phenotype_and_metric_from_file(n, dtype=float, ignore="............")
        ph_p, peaks_per_fl =  load_phenotype_and_metric_from_file(p, dtype=int, ignore="............")

        y_values = []
        peak_ratios_per_fl = []
        local_peak_frac_per_fl = []
        # loop over fitness landscapes
        for glo_peak_ph, navig, all_peaks in zip(ph_n, navigs_per_fl, peaks_per_fl):
            peak_ratio = ph_to_counts[glo_peak_ph] / all_peaks  # global peaks / all peaks
            peak_ratios_per_fl.append(peak_ratio)
            local_peak_frac_per_fl.append((all_peaks - ph_to_counts[glo_peak_ph]) / folded_count)
            
        ax1.scatter(peak_ratios_per_fl, navigs_per_fl, color=f"C{i-1}", alpha=0.5, s=7, linewidths=0)

        mean_peak_ratio = np.mean(peak_ratios_per_fl)
        mean_navig = np.mean(navigs_per_fl)

        mean_peak_ratios.append(mean_peak_ratio)
        mean_navigs.append(mean_navig)

        xq1 = [np.max([mean_peak_ratio-np.percentile(peak_ratios_per_fl, q=25), 0])]
        xq2 = [np.max([np.percentile(peak_ratios_per_fl, q=75)-mean_peak_ratio, 0])]
        
        yq1 = [np.max([mean_navig-np.percentile(navigs_per_fl, q=25), 0])]
        yq2 = [np.max([np.percentile(navigs_per_fl, q=75)-mean_navig, 0])]

        ax1.errorbar(mean_peak_ratio, mean_navig, xerr=(xq1, xq2), yerr=(yq1, yq2), elinewidth=1, marker="s", markersize=5, label=f"GP map {i}", color=f"C{i-1}")
    
        ax1.plot([0, 1], [0, 1], color="grey", linestyle="dotted", linewidth=0.5, zorder=-10)

        mean_loc_peak_frac = np.mean(local_peak_frac_per_fl)
      
        mean_loc_peak_fracs.append(mean_loc_peak_frac)
    
        xq1 = [np.max([mean_loc_peak_frac-np.percentile(local_peak_frac_per_fl, q=25), 0])]
        xq2 = [np.max([np.percentile(local_peak_frac_per_fl, q=75)-mean_loc_peak_frac, 0])]

        ax2.errorbar(mean_loc_peak_frac, mean_navig, xerr=(xq1, xq2), yerr=(yq1, yq2), elinewidth=1, marker="s", markersize=5, label=f"GP map {i}", color=f"C{i-1}")

    r, p = pearsonr(mean_peak_ratios, mean_navigs)
    p_str = "%.3g" % p
    ax1.text(.04, .9, f'r = {np.round(r, 2)}\np = {p_str}', transform=ax1.transAxes, horizontalalignment='left', size=12)
    
    ax1.set_xlabel("f_glo/(f_loc+f_glo)", size=15)
    ax1.set_ylabel("Navigability", size=15)
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1)
    ax1.tick_params(axis='both', which='major', labelsize=15)
    ax1.tick_params(axis='both', which='minor', labelsize=15)
    # ax1.legend(fontsize=6, title="GP map", loc="lower right", frameon=False)

    r, p = pearsonr(mean_loc_peak_fracs, mean_navigs)
    p_str = "%.3g" % p
    ax2.text(.6, .9, f'r = {np.round(r, 2)}\np = {p_str}', transform=ax2.transAxes, horizontalalignment='left', size=12)
    
    ax2.set_xlabel("local peak fraction f_loc", size=15)
    ax2.set_ylabel("Navigability", size=15)
    ax2.set_ylim(0, 1)
    ax2.tick_params(axis='both', which='major', labelsize=15)
    ax2.tick_params(axis='both', which='minor', labelsize=15)

    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)
