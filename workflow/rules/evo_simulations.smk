

rule build_random_fitness_landscape:
    input:
        "<output>/gp_map_{bp}/phenotypes.txt"
    output:
        "<output>/gp_map_{bp}/fitness_landscapes/fl_{fl_id}/fl.txt"
    params:
        min_f=config["fl_params"]["min_f"],
        max_f=config["fl_params"]["max_f"],
        lethal_ph=config["unfolded"]
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    shell:
        "workflow/scripts/evo_simulations/random_fitness_landscape.py "
        "--phenotypes {input} "
        "--low_f {params.min_f} "
        "--upp_f {params.max_f} "
        "--lethal_ph {params.lethal_ph} "
        "-o {output} "

rule run_adaptive_walks_productive_trace_paths:
    input:
        fl="<output>/gp_map_{bp}/fitness_landscapes/fl_{fl_id}/fl.txt",
        gp_map="<output>/gp_map_{bp}/gp_map.pickle"
    output:
        walk_lengths="<output>/gp_map_{bp}/fitness_landscapes/fl_{fl_id}/adaptive_walks/pop_size{pop_size}/walk_lengths_seed.csv",
        paths="<output>/gp_map_{bp}/fitness_landscapes/fl_{fl_id}/adaptive_walks/pop_size{pop_size}/paths_seed.csv"
    params:        
        max_steps=config["evo_sim_params"]["max_steps"],  
        sample_size_walks=config["evo_sim_params"]["sample_size"],
        avoid_phenotype=config["unfolded"]
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    shell:
        "workflow/scripts/evo_simulations/adaptive_walks_productive_trace_path.py "
        "-i {input.gp_map} "
        "-f {input.fl} "
        "-l {output.walk_lengths} "
        "-p {output.paths} "
        "-n {wildcards.pop_size} "
        "-m {params.max_steps} "
        "-a {params.avoid_phenotype} "
        "-s {params.sample_size_walks} "

rule compute_navigability_per_fl:
    input:
        walk_length_files=expand("<output>/gp_map_{{bp}}/fitness_landscapes/fl_{fl_id}/adaptive_walks/pop_size{{pop_size}}/walk_lengths_seed.csv", fl_id=range(config["fl_params"]["N"])),
        fitness_landscapes=expand("<output>/gp_map_{{bp}}/fitness_landscapes/fl_{fl_id}/fl.txt", fl_id=range(config["fl_params"]["N"])),
    output:
        "<output>/gp_map_{bp}/navigability_per_fl_pop{pop_size}.txt"
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    shell:
        "workflow/scripts/evo_simulations/compute_navigability_per_fl.py "
        "-w {input.walk_length_files} "  # assumes one per fitness landscape
        "-f {input.fitness_landscapes} "
        "-o {output} "
