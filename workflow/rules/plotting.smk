
rule plot_boxplots_navigability:
    input:
        expand("<output>/gp_map_{bp}/navigability_per_fl_pop{{pop_size}}.txt", bp=config["bp_rules"])
    output:
        "<output>/figures/navigability_boxplot_pop{pop_size}.pdf"
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    shell:
        "workflow/scripts/plotting/plot_navigability_boxplots.py "
        "-n {input} "
        "-o {output} "

rule plot_vienna_comparison_triple_plot:
    input:
        ph_distr="<output>/gp_map_{bp}/phenotype_distribution.txt",
        ref_ph_distr="<output>/gp_map_{ref_map}/phenotype_distribution.txt",
        nc="<output>/gp_map_{bp}/neutral_component_sizes.txt",
        ref_nc="<output>/gp_map_{ref_map}/neutral_component_sizes.txt",
        walks="<output>/gp_map_{bp}/navigability_per_fl_pop{pop_size}.txt",
        ref_walks="<output>/gp_map_{ref_map}/navigability_per_fl_pop{pop_size}.txt",
    output:
        "<output>/figures/{ref_map}_bp{bp}_comparison_pop{pop_size}.pdf"
    params:
        sample_size=config["evo_sim_params"]["sample_size"],  # how many walks per phenotype per random fitness landscape instance
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    shell:
        "workflow/scripts/plotting/plot_viennaRNA_comparison_triple_plot.py "
        "-i {input.ph_distr} "
        "-r {input.ref_ph_distr} "
        "-n {input.nc} "
        "-m {input.ref_nc} "
        "--refwalks {input.ref_walks} "
        "--walks {input.walks} "
        "--sample_size {params.sample_size} "
        "-o {output} "

rule plot_triple_comp_navig_over_peak_freq:
    input:
        ph_dists=["<output>/gp_map_{canon}/phenotype_distribution.txt", "<output>/gp_map_{comp1}/phenotype_distribution.txt", "<output>/gp_map_{comp2}/phenotype_distribution.txt"],
        navigs=["<output>/gp_map_{canon}/navigability_per_fl_pop{pop_size}.txt", "<output>/gp_map_{comp1}/navigability_per_fl_pop{pop_size}.txt", "<output>/gp_map_{comp2}/navigability_per_fl_pop{pop_size}.txt"]
    output:
        "<output>/figures/navigability_over_freq_bp{canon}_bp{comp1}_bp{comp2}_pop_size{pop_size}.pdf"
    params:
        ignore=config["unfolded"]
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    shell:
        "workflow/scripts/plotting/plot_triple_comp_navig_over_peak_freq.py "
        "-p {input.ph_dists} "
        "-n {input.navigs} "
        "-i {params.ignore} "
        "-o {output} "

rule plot_navig_over_peak_freq_all:
    input:
        ph_dists=expand("<output>/gp_map_{bp}/phenotype_distribution.txt", bp=config["bp_rules"]),
        navigs=expand("<output>/gp_map_{bp}/navigability_per_fl_pop{{pop_size}}.txt", bp=config["bp_rules"])
    output:
        "<output>/figures/navigability_over_freq_all_pop_size{pop_size}.pdf"
    params:
        ignore=config["unfolded"]
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    shell:
        "workflow/scripts/plotting/plot_navig_over_peak_freq_all.py "
        "-p {input.ph_dists} "
        "-n {input.navigs} "
        "-i {params.ignore} "
        "-o {output} "

rule plot_navigability_over_peak_ratio:
    input:
        per_fl_navigability=expand("<output>/gp_map_{bp}/navigability_per_fl_pop{{pop_size}}.txt", bp=config["bp_rules"]),
        phenotype_distribution=expand("<output>/gp_map_{bp}/phenotype_distribution.txt", bp=config["bp_rules"]),
        per_fl_combined_peak_size=expand("<output>/gp_map_{bp}/combined_peak_sizes_per_fl.txt", bp=config["bp_rules"]),
    output:
        "<output>/figures/navigability_over_peak_ratio_pop{pop_size}.pdf"
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    shell:
        "workflow/scripts/plotting/plot_navigability_over_peak_ratio.py "
        "-n {input.per_fl_navigability} "
        "-f {input.phenotype_distribution} "
        "-p {input.per_fl_combined_peak_size} "
        "-o {output} "

rule plot_peak_sizes_vs_prediction:
    input:
        nc_graphs=expand("<output>/gp_map_{bp}/nc_graph.pickle", bp=config["bp_rules"]),
        per_fl_combined_peak_size=expand("<output>/gp_map_{bp}/combined_peak_sizes_per_fl.txt", bp=config["bp_rules"]),
    output:
        "<output>/figures/peak_sizes_vs_prediction{pop_size}.pdf"
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    shell:
        "workflow/scripts/plotting/plot_peak_sizes_vs_prediction.py "
        "-n {input.nc_graphs} "
        "-p {input.per_fl_combined_peak_size} "
        "-o {output} "

rule plot_neigh_div_degeneracy_relation:
    input:
        nc_to_gt_files=expand("<output>/gp_map_{bp}/nc_to_gt.txt", bp=config["bp_rules"]),
        gp_maps=expand("<output>/gp_map_{bp}/gp_map.pickle", bp=config["bp_rules"])
    output:
        "<output>/figures/fig6_neigh_div_degen_relation.pdf"
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    shell:
        "workflow/scripts/plotting/plot_neigh_div_bp_deg_relation.py "
        "-g {input.gp_maps} "
        "-n {input.nc_to_gt_files} "
        "-o {output} "