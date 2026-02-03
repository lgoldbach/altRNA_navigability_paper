
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
