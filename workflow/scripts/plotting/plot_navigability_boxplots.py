#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
from matplotlib import colormaps
from matplotlib.ticker import MaxNLocator
import numpy as np
import seaborn as sns

from rna_gpf.parsing import read_navigability_per_fl


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--navigability", help="Input files with " \
                        "navigability values for each phenotype for each fitness landscape"
                        "format: (((...))) 5 2 1 4 2 4 ", 
                        required=True, nargs='+')
    parser.add_argument("-o", "--output", help="Output file name "
                        "(should end in .pdf)", required=True)
    
    args = parser.parse_args()
    fig, ax = plt.subplots()

    
    navigs = [read_navigability_per_fl(file) for file in args.navigability]

    # get navigabilities in seaborn usable format
    navigs_per_bp = []
    navigs_sns = []
    x = []
    for i, nav in enumerate(navigs, start=1):
        navigs_per_bp.append([])
        for navigs_over_fls in nav.values():  # per ph navigs
            for n in navigs_over_fls:
                x.append(i)
                navigs_sns.append(n)
                navigs_per_bp[-1].append(n)
            

    palette = ["0.9" if i != 4 else "0.6" for i in range(1, len(x))]
    ax = sns.boxplot(x=x,
                y=navigs_sns,
                ax=ax,
                color="0.9",
                linewidth=0.2,
                linecolor="black",
                legend=False,
                showfliers=False,
                zorder=3,
                palette=palette,
                whis=[5, 95])  

    # ax = sns.violinplot(x=x,
    #             y=navigs_sns,
    #             palette=palette,
    #             ax=ax,
    #             legend=False,
    #             density_norm="area",
    #             width=0.95,
    #             common_norm=True,
    #             cut=0,
    #             inner="quart",
    #             linewidth=0.2,
    #             linecolor="black",
    #             zorder=3,
    #             inner_kws={"zorder": 4}) 
    
    for l in ax.lines:
        l.set_linestyle('-')
        l.set_color('black')
        l.set_linewidth(0.5)


    for i, navi in enumerate(navigs_per_bp):
        d = navi
        q1 = np.percentile(d, q=25)
        q3 = np.percentile(d, q=75)
        mean = np.mean(d)
        median = np.median(d)

        print(i, mean)
        # ax1.vlines(i, q1, q3, color="black", linewidth=1, zorder=5)
        if i == 1:  # only add label once for legend
            ax.plot([i-.1, i+0.1], [mean, mean], color="red", zorder=10, linewidth=3, label="Mean: g-p map navigability") 
        else:
            ax.plot([i-.1, i+0.1], [mean, mean], color="red", zorder=10, linewidth=3)


    ax.set_ylabel("Fraction of successful adaptive walks\nper target phenotype", fontsize=15)
    ax.set_xlabel("g-p map", fontsize=15)

    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.tick_params(axis='both', which='minor', labelsize=15)

    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xticks(list(range(0, 10)))
    p =[1, 2, 3, "canon.", 5, 6, 7, 8, 9, 10]
    ax.set_xticklabels(p) #(list(range(1, 11)))

    ax.grid(axis="y", zorder=-1)

    # l = ax2.legend(title="Selection pressure", prop={'size': 12})
    # plt.setp(l.get_title(),fontsize=12)

    # ax.legend(loc="lower left", frameon=False, fancybox=False, prop={'size': 15})

    ax.set_ylim(0, 1)

    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)
