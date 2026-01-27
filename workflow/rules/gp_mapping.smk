

rule generate_alt_genotype_space:
    output:
        "<output>/alt_genotypes.txt"
    params:
        alphabet=config["alt_alphabet"],
        seq_len=config["sequence_length"]
    shell:
        "workflow/scripts/gp_mapping/build_genotype_space.py " 
        "-a {params.alphabet} "
        "-l {params.seq_len} "
        "-o {output}"

rule generate_viennaRNA_genotype_space:
    output:
        "<output>/gp_map_viennaRNA/genotypes.txt"
    params:
        alphabet=config["viennaRNA"]["alphabet"],
        seq_len=config["sequence_length"]
    shell:
        "workflow/scripts/gp_mapping/build_genotype_space.py " 
        "-a {params.alphabet} "
        "-l {params.seq_len} "
        "-o {output}"

rule make_permissible_sets:
    input:
        "<output>/alt_genotypes.txt"
    output:
        "workflow/scripts/gp_map_{bp}/permissible_sets.txt"  # regex to constrain to numbers and prevent matching subdirectory gp_map.txt
    params:
        graph_path=config["bp_rule_graphs"],
        min_loop_size=config["folding_params"]["min_loop"],
        suboptimal=config["folding_params"]["hbonds_below_max"],
        structures_max=config["folding_params"]["permis_set_max"],
        alphabet=config["alt_alphabet"]
    shell:
        "workflow/scripts/gp_mapping/generate_permissible_sets.py "
        "-i {input} "
        "-o {output} "
        "-m {params.min_loop_size} "
        "-s {params.suboptimal} "
        "-z {params.structures_max} "
        "-g {params.graph_path} "
        "-a {params.alphabet} "
        "-p {wildcards.bp} "

rule viennaRNA_mfe_gp_map:
    input:
        "<output>/gp_map_viennaRNA/genotypes.txt"
    output:
        "<output>/gp_map_viennaRNA/gp_map.txt"
    shell:
        "workflow/scripts/gp_mapping/viennaRNAmfe_gp_mapping.py "
        "-i {input} "
        "-o {output} "