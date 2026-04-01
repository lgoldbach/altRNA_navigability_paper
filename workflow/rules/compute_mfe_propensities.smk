
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
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    shell:
        "workflow/scripts/compute_mfe_propensities/preprocess_gp_map_and_permissible_sets.py "
        "-i {input.permissible_sets} "
        "-g1 {input.permissible_set_genotypes} "
        "-r {input.gp_map} "
        "-g2 {input.gp_map_genotypes} "
        "-d {params.dead_ph} "
        "-o1 {output.output_subopt} "
        "-o2 {output.output_mfe} "

rule pairwise_consensus_matrix:
    input:
        gp_subopt="<output>/gp_map_" + config["mfe_prop_ref_gp_map"] + "/mfe_prop/pg_map_supoptimal_augc_folded_dict.pickle",
        gp_ref="<output>/gp_map_" + config["mfe_prop_ref_gp_map"] + "/mfe_prop/pg_map_mfe_augc_folded_dict.pickle"
    output:
        matrix="<output>/gp_map_" + config["mfe_prop_ref_gp_map"] + "/pairwise_consensus_matrix.pickle",
        phenotype_list="<output>/gp_map_" + config["mfe_prop_ref_gp_map"] + "/mfe_ranking_phenotypes.txt"
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    shell:
        "workflow/scripts/compute_mfe_propensities/pairwise_consensus_matrix.py "
        "-i {input.gp_subopt} "
        "-r {input.gp_ref} "
        "-o {output.matrix} "
        "-p {output.phenotype_list} "

rule compute_mfe_propensities:
    input:
        matrix="<output>/gp_map_" + config["mfe_prop_ref_gp_map"] + "/pairwise_consensus_matrix.pickle",
        phenotypes="<output>/gp_map_" + config["mfe_prop_ref_gp_map"] + "/mfe_ranking_phenotypes.txt"
    output:
        "<output>/gp_map_" + config["mfe_prop_ref_gp_map"] + "/mfe_propensities.csv"
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    shell:
        "workflow/scripts/compute_mfe_propensities/bradley_terry_ranking.py "
        "-i {input.matrix} "
        "-p {input.phenotypes} "
        "-o {output} "


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
