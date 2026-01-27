
script_dir = ""

rule preprocess_suboptimal_map:
    input: 
        permissible_sets="<output>/gp_map_" + config["canonical_rule"] + "/permissible_sets.txt",
        permissible_set_genotypes="<output>/alt_genotypes.txt",
        gp_map="<output>/gp_map_" + config["mfe_prop_ref_gp_map"] + "/gp_map.txt",
        gp_map_genotypes="<output>/gp_map_" + config["mfe_prop_ref_gp_map"] + "/genotypes.txt"
    output: 
        output_subopt="<output>/gp_map_" + config["mfe_prop_ref_gp_map"] + "/mfe_prop/pg_map_supoptimal_augc_folded_dict.pickle",
        output_mfe="<output>/gp_map_" + config["mfe_prop_ref_gp_map"] + "/mfe_prop/pg_map_mfe_augc_folded_dict.pickle"
    params:
        dead_ph=config["unfolded"]
    shell:
        "workflow/scripts//mfe_propensity_analysis/prepocess_suboptiomal_map.py "
        "-i {input.permissible_sets} "
        "-g1 {input.permissible_set_genotypes} "
        "-r {input.gp_map} "
        "-g2 {input.gp_map_genotypes} "
        "-d {params.dead_ph} "
        "-o1 {output.output_subopt} "
        "-o2 {output.output_mfe} "

# rule pairwise_consensus_matrix:
#     input:
#         gp_subopt=config["ref_gp_map_suboptimal"],
#         gp_ref=config["ref_gp_map"]
#     output:
#         matrix="pairwise_consensus_matrix.pickle",
#         phenotype_list="phenotypes.txt"
#     shell:
#         "pairwise_consensus_matrix.py "
#         "-i {input.pg_subopt} "
#         "-r {input.gp_ref} "
#         "-o {output.matrix} "
#         "-p {output.phenotype_list} "
    
# rule consenus_matrix_statistics:
#     input:
#         "pairwise_consensus_matrix.pickle"
#     output:
#         "consistency_matrix_stats.pdf"
#     shell:
#         "plot_consensus_matrix_stats.py "
#         "-i {input} "
#         "-o {output} "

# rule binarize_consensus_matrix:
#     input:
#         "pairwise_consensus_matrix.pickle"
#     output:
#         "binary_pairwise_consensus_matrix.pickle"
#     shell:
#         "binarize_consensus_matrix.py "
#         "-i {input} "
#         "-o {output} "

# rule plot_consensus_matrix:
#     input:
#         matrix="pairwise_consensus_matrix.pickle",
#         phenotypes="phenotypes.txt"
#     output:
#         "consensus_matrix_plot.pdf"
#     shell:
#         "plot_consensus_matrix.py "
#         "-i {input.matrix} "
#         "-p {input.phenotypes} "
#         "-o {output} "

# rule topological_sorting:
#     input:
#         matrix="binary_pairwise_consensus_matrix.pickle",
#         phenotypes="phenotypes.txt"
#     output:
#         "rankings.csv"
#     shell:
#         "topological_sorting.py "
#         "-i {input.matrix} "
#         "-o {output} "
#         "-p {input.phenotypes} "

    

