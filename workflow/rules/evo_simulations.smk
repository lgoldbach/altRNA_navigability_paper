

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
    wildcard_constraints:
        pop_size='[0-9]*'
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
    wildcard_constraints:
         pop_size='[0-9]*'
    shell:
        "workflow/scripts/evo_simulations/compute_navigability_per_fl.py "
        "-w {input.walk_length_files} "  # assumes one per fitness landscape
        "-f {input.fitness_landscapes} "
        "-o {output} "

rule compute_combined_peak_size_per_fl:
    input:
        nc_graph="<output>/gp_map_{bp}/nc_graph.pickle",
        fitness_landscapes=expand("<output>/gp_map_{{bp}}/fitness_landscapes/fl_{fl_id}/fl.txt", fl_id=range(config["fl_params"]["N"])),
    output:
        "<output>/gp_map_{bp}/combined_peak_sizes_per_fl.txt"
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    shell:
        "workflow/scripts/evo_simulations/compute_combined_peak_size_per_fl.py "
        "-n {input.nc_graph} "
        "-f {input.fitness_landscapes} "
        "-o {output} "



#### correlated landscapes

rule build_correlated_fitness_landscape:
    input:
        "<output>/gp_map_{bp}/phenotypes.txt"
    output:
        "<output>/gp_map_{bp}/correlated_fitness_landscapes/cparam_{cparam}/fl_{fl_id}/fl.txt"
    params:
        min_f=config["fl_params"]["min_f"],
        max_f=config["fl_params"]["max_f"],
        lethal_ph=config["unfolded"]
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    shell:
        "workflow/scripts/evo_simulations/correlated_fitness_landscape.py "
        "--phenotype_file {input} "
        "--decay linear "
        "--mix {wildcards.cparam} "
        "--output {output} "

rule run_adaptive_walks_productive_trace_paths_corr:
    input:
        fl="<output>/gp_map_{bp}/correlated_fitness_landscapes/cparam_{cparam}/fl_{fl_id}/fl.txt",
        gp_map="<output>/gp_map_{bp}/gp_map.pickle"
    output:
        walk_lengths="<output>/gp_map_{bp}/correlated_fitness_landscapes/cparam_{cparam}/fl_{fl_id}/adaptive_walks/pop_size{pop_size}/walk_lengths_seed.csv",
        paths="<output>/gp_map_{bp}/correlated_fitness_landscapes/cparam_{cparam}/fl_{fl_id}/adaptive_walks/pop_size{pop_size}/paths_seed.csv"
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

rule compute_navigability_per_fl_corr:
    input:
        walk_length_files=expand("<output>/gp_map_{{bp}}/correlated_fitness_landscapes/cparam_{{cparam}}/fl_{fl_id}/adaptive_walks/pop_size{{pop_size}}/walk_lengths_seed.csv", fl_id=range(config["corr_fl_params"]["N"])),
        fitness_landscapes=expand("<output>/gp_map_{{bp}}/correlated_fitness_landscapes/cparam_{{cparam}}/fl_{fl_id}/fl.txt", fl_id=range(config["corr_fl_params"]["N"])),
    output:
        "<output>/gp_map_{bp}/navigability_per_fl_pop{pop_size}_corr{cparam}.txt"
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    shell:
        "workflow/scripts/evo_simulations/compute_navigability_per_fl.py "
        "-w {input.walk_length_files} "  # assumes one per fitness landscape
        "-f {input.fitness_landscapes} "
        "-o {output} "


### random unfolded
rule build_random_fitness_landscape_random_unfolded:
    input:
        "<output>/gp_map_{bp}/phenotypes_unfolded{n_unfolded}.txt"
    output:
        "<output>/gp_map_{bp}/unfolded{n_unfolded}/fitness_landscapes/fl_{fl_id}/fl.txt"
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

rule run_adaptive_walks_productive_trace_paths_random_unfolded:
    input:
        fl="<output>/gp_map_{bp}/unfolded{n_unfolded}/fitness_landscapes/fl_{fl_id}/fl.txt",
        gp_map="<output>/gp_map_{bp}/gp_map_unfolded{n_unfolded}.pickle"
    output:
        walk_lengths="<output>/gp_map_{bp}/unfolded{n_unfolded}/fitness_landscapes/fl_{fl_id}/adaptive_walks/pop_size{pop_size}/walk_lengths_seed.csv",
        paths="<output>/gp_map_{bp}/unfolded{n_unfolded}/fitness_landscapes/fl_{fl_id}/adaptive_walks/pop_size{pop_size}/paths_seed.csv"
    params:        
        max_steps=config["evo_sim_params"]["max_steps"],  
        sample_size_walks=config["evo_sim_params"]["sample_size"],
        avoid_phenotype=config["unfolded"]
    wildcard_constraints:
        pop_size='[0-9]*'
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

rule compute_navigability_per_fl_random_unfolded:
    input:
        walk_length_files=expand("<output>/gp_map_{{bp}}/unfolded{{n_unfolded}}/fitness_landscapes/fl_{fl_id}/adaptive_walks/pop_size{{pop_size}}/walk_lengths_seed.csv", fl_id=range(config["fl_params"]["N"])),
        fitness_landscapes=expand("<output>/gp_map_{{bp}}/unfolded{{n_unfolded}}/fitness_landscapes/fl_{fl_id}/fl.txt", fl_id=range(config["fl_params"]["N"])),
    output:
        "<output>/gp_map_{bp}/navigability_per_fl_pop{pop_size}_unfolded{n_unfolded}.txt"
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    wildcard_constraints:
         pop_size='[0-9]*'
    shell:
        "workflow/scripts/evo_simulations/compute_navigability_per_fl.py "
        "-w {input.walk_length_files} "  # assumes one per fitness landscape
        "-f {input.fitness_landscapes} "
        "-o {output} "


### ref unfolded
rule build_random_fitness_landscape_ref_unfolded:
    input:
        "<output>/gp_map_{bp}/phenotypes_ref_unfolded.txt"
    output:
        "<output>/gp_map_{bp}/ref_unfolded/fitness_landscapes/fl_{fl_id}/fl.txt"
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

rule run_adaptive_walks_productive_trace_paths_ref_unfolded:
    input:
        fl="<output>/gp_map_{bp}/ref_unfolded/fitness_landscapes/fl_{fl_id}/fl.txt",
        gp_map="<output>/gp_map_{bp}/gp_map_ref_unfolded.pickle"
    output:
        walk_lengths="<output>/gp_map_{bp}/ref_unfolded/fitness_landscapes/fl_{fl_id}/adaptive_walks/pop_size{pop_size}/walk_lengths_seed.csv",
        paths="<output>/gp_map_{bp}/ref_unfolded/fitness_landscapes/fl_{fl_id}/adaptive_walks/pop_size{pop_size}/paths_seed.csv"
    params:        
        max_steps=config["evo_sim_params"]["max_steps"],  
        sample_size_walks=config["evo_sim_params"]["sample_size"],
        avoid_phenotype=config["unfolded"]
    wildcard_constraints:
        pop_size='[0-9]*'
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

rule compute_navigability_per_fl_ref_unfolded:
    input:
        walk_length_files=expand("<output>/gp_map_{{bp}}/ref_unfolded/fitness_landscapes/fl_{fl_id}/adaptive_walks/pop_size{{pop_size}}/walk_lengths_seed.csv", fl_id=range(config["fl_params"]["N"])),
        fitness_landscapes=expand("<output>/gp_map_{{bp}}/ref_unfolded/fitness_landscapes/fl_{fl_id}/fl.txt", fl_id=range(config["fl_params"]["N"])),
    output:
        "<output>/gp_map_{bp}/navigability_per_fl_pop{pop_size}_ref_unfolded.txt"
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    wildcard_constraints:
         pop_size='[0-9]*'
    shell:
        "workflow/scripts/evo_simulations/compute_navigability_per_fl.py "
        "-w {input.walk_length_files} "  # assumes one per fitness landscape
        "-f {input.fitness_landscapes} "
        "-o {output} "


# nc unfolded
rule build_random_fitness_landscape_nc_unfolded:
    input:
        "<output>/gp_map_{bp}/phenotypes_nc_unfolded_n_unf{n_unfolded}.txt",
    output:
        "<output>/gp_map_{bp}/nc_unfolded_{n_unfolded}/fitness_landscapes/fl_{fl_id}/fl.txt"
    params:
        min_f=config["fl_params"]["min_f"],
        max_f=config["fl_params"]["max_f"],
        lethal_ph=config["unfolded"]
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    wildcard_constraints:
        bp='[0-9]*'
    shell:
        "workflow/scripts/evo_simulations/random_fitness_landscape.py "
        "--phenotypes {input} "
        "--low_f {params.min_f} "
        "--upp_f {params.max_f} "
        "--lethal_ph {params.lethal_ph} "
        "-o {output} "


rule run_adaptive_walks_productive_trace_paths_nc_unfolded:
    input:
        fl="<output>/gp_map_{bp}/nc_unfolded_{n_unfolded}/fitness_landscapes/fl_{fl_id}/fl.txt",
        gp_map="<output>/gp_map_{bp}/gp_map_nc_unfolded_n_unf{n_unfolded}.pickle",
    output:
        walk_lengths="<output>/gp_map_{bp}/nc_unfolded_{n_unfolded}/fitness_landscapes/fl_{fl_id}/adaptive_walks/pop_size{pop_size}/walk_lengths_seed.csv",
        paths="<output>/gp_map_{bp}/nc_unfolded_{n_unfolded}/fitness_landscapes/fl_{fl_id}/adaptive_walks/pop_size{pop_size}/paths_seed.csv"
    params:        
        max_steps=config["evo_sim_params"]["max_steps"],  
        sample_size_walks=config["evo_sim_params"]["sample_size"],
        avoid_phenotype=config["unfolded"]
    wildcard_constraints:
        pop_size='[0-9]*',
        bp='[0-9]*'
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

rule compute_navigability_per_fl_nc_unfolded:
    input:
        walk_length_files=expand("<output>/gp_map_{{bp}}/nc_unfolded_{{n_unfolded}}/fitness_landscapes/fl_{fl_id}/adaptive_walks/pop_size{{pop_size}}/walk_lengths_seed.csv", fl_id=range(config["fl_params"]["N"])),
        fitness_landscapes=expand("<output>/gp_map_{{bp}}/nc_unfolded_{{n_unfolded}}/fitness_landscapes/fl_{fl_id}/fl.txt", fl_id=range(config["fl_params"]["N"])),
    output:
        "<output>/gp_map_{bp}/navigability_per_fl_pop{pop_size}_nc_unfolded_{n_unfolded}.txt"
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    wildcard_constraints:
        pop_size='[0-9]*',
        bp='[0-9]*'
    shell:
        "workflow/scripts/evo_simulations/compute_navigability_per_fl.py "
        "-w {input.walk_length_files} "  # assumes one per fitness landscape
        "-f {input.fitness_landscapes} "
        "-o {output} "




# diffu unfolded
rule build_random_fitness_landscape_diffu_unfolded:
    input:
        "<output>/gp_map_{bp}/phenotypes_diffu_unfolded.txt"
    output:
        "<output>/gp_map_{bp}/diffu/fitness_landscapes/fl_{fl_id}/fl.txt"
    params:
        min_f=config["fl_params"]["min_f"],
        max_f=config["fl_params"]["max_f"],
        lethal_ph=config["unfolded"]
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    wildcard_constraints:
        bp='[0-9]*'
    shell:
        "workflow/scripts/evo_simulations/random_fitness_landscape.py "
        "--phenotypes {input} "
        "--low_f {params.min_f} "
        "--upp_f {params.max_f} "
        "--lethal_ph {params.lethal_ph} "
        "-o {output} "


rule run_adaptive_walks_productive_trace_paths_diffu_unfolded:
    input:
        fl="<output>/gp_map_{bp}/diffu/fitness_landscapes/fl_{fl_id}/fl.txt",
        gp_map="<output>/gp_map_{bp}/gp_map_diffu_unfolded.pickle",
    output:
        walk_lengths="<output>/gp_map_{bp}/diffu/fitness_landscapes/fl_{fl_id}/adaptive_walks/pop_size{pop_size}/walk_lengths_seed.csv",
        paths="<output>/gp_map_{bp}/diffu/fitness_landscapes/fl_{fl_id}/adaptive_walks/pop_size{pop_size}/paths_seed.csv"
    params:        
        max_steps=config["evo_sim_params"]["max_steps"],  
        sample_size_walks=config["evo_sim_params"]["sample_size"],
        avoid_phenotype=config["unfolded"]
    wildcard_constraints:
        pop_size='[0-9]*',
        bp='[0-9]*'
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

rule compute_navigability_per_fl_diffu_unfolded:
    input:
        walk_length_files=expand("<output>/gp_map_{{bp}}/diffu/fitness_landscapes/fl_{fl_id}/adaptive_walks/pop_size{{pop_size}}/walk_lengths_seed.csv", fl_id=range(config["fl_params"]["N"])),
        fitness_landscapes=expand("<output>/gp_map_{{bp}}/diffu/fitness_landscapes/fl_{fl_id}/fl.txt", fl_id=range(config["fl_params"]["N"])),
    output:
        "<output>/gp_map_{bp}/navigability_per_fl_pop{pop_size}_diffu_unfolded.txt"
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    wildcard_constraints:
        pop_size='[0-9]*',
        bp='[0-9]*'
    shell:
        "workflow/scripts/evo_simulations/compute_navigability_per_fl.py "
        "-w {input.walk_length_files} "  # assumes one per fitness landscape
        "-f {input.fitness_landscapes} "
        "-o {output} "



# ph based unfolded
rule build_random_fitness_landscape_ph_based_unfolded:
    input:
        "<output>/gp_map_{bp}/phenotypes_ph_based_unfolded.txt"
    output:
        "<output>/gp_map_{bp}/ph_based_unfolded/fitness_landscapes/fl_{fl_id}/fl.txt"
    params:
        min_f=config["fl_params"]["min_f"],
        max_f=config["fl_params"]["max_f"],
        lethal_ph=config["unfolded"]
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    wildcard_constraints:
        bp='[0-9]*'
    shell:
        "workflow/scripts/evo_simulations/random_fitness_landscape.py "
        "--phenotypes {input} "
        "--low_f {params.min_f} "
        "--upp_f {params.max_f} "
        "--lethal_ph {params.lethal_ph} "
        "-o {output} "

rule run_adaptive_walks_productive_trace_paths_ph_based_unfolded:
    input:
        fl="<output>/gp_map_{bp}/ph_based_unfolded/fitness_landscapes/fl_{fl_id}/fl.txt",
        gp_map="<output>/gp_map_{bp}/gp_map_ph_based_unfolded.pickle",
    output:
        walk_lengths="<output>/gp_map_{bp}/ph_based_unfolded/fitness_landscapes/fl_{fl_id}/adaptive_walks/pop_size{pop_size}/walk_lengths_seed.csv",
        paths="<output>/gp_map_{bp}/ph_based_unfolded/fitness_landscapes/fl_{fl_id}/adaptive_walks/pop_size{pop_size}/paths_seed.csv"
    params:        
        max_steps=config["evo_sim_params"]["max_steps"],  
        sample_size_walks=config["evo_sim_params"]["sample_size"],
        avoid_phenotype=config["unfolded"]
    wildcard_constraints:
        pop_size='[0-9]*',
        bp='[0-9]*'
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

rule compute_navigability_per_fl_ph_based_unfolded:
    input:
        walk_length_files=expand("<output>/gp_map_{{bp}}/ph_based_unfolded/fitness_landscapes/fl_{fl_id}/adaptive_walks/pop_size{{pop_size}}/walk_lengths_seed.csv", fl_id=range(config["fl_params"]["N"])),
        fitness_landscapes=expand("<output>/gp_map_{{bp}}/ph_based_unfolded/fitness_landscapes/fl_{fl_id}/fl.txt", fl_id=range(config["fl_params"]["N"])),
    output:
        "<output>/gp_map_{bp}/navigability_per_fl_pop{pop_size}_ph_based_unfolded.txt"
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    wildcard_constraints:
        pop_size='[0-9]*',
        bp='[0-9]*'
    shell:
        "workflow/scripts/evo_simulations/compute_navigability_per_fl.py "
        "-w {input.walk_length_files} "  # assumes one per fitness landscape
        "-f {input.fitness_landscapes} "
        "-o {output} "


# comparing navigability seq-free and viennaRNA
rule run_adaptive_walks_productive_trace_paths_seq_same_as_vienna:
    input:
        fl="<output>/gp_map_04/fitness_landscapes/fl_{fl_id}/fl.txt",  # use gp_map 4 fitness landscapes 
        gp_map="<output>/gp_map_viennaRNA/gp_map.pickle",  # and viennaRNA gp map
    output:
        walk_lengths="<output>/gp_map_04/vienna_RNA_gpmap/fitness_landscapes/fl_{fl_id}/adaptive_walks/pop_size{pop_size}/walk_lengths_seed.csv",
        paths="<output>/gp_map_04/vienna_RNA_gpmap/fitness_landscapes/fl_{fl_id}/adaptive_walks/pop_size{pop_size}/paths_seed.csv"
    params:        
        max_steps=config["evo_sim_params"]["max_steps"],  
        sample_size_walks=config["evo_sim_params"]["sample_size"],
        avoid_phenotype=config["unfolded"]
    wildcard_constraints:
        pop_size='[0-9]*',
        bp='[0-9]*'
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

rule compute_navigability_per_fl_seq_same_as_vienna:
    input:
        walk_length_files=expand("<output>/gp_map_04/vienna_RNA_gpmap/fitness_landscapes/fl_{fl_id}/adaptive_walks/pop_size{{pop_size}}/walk_lengths_seed.csv", fl_id=range(config["fl_params"]["N"])),
        fitness_landscapes=expand("<output>/gp_map_04/fitness_landscapes/fl_{fl_id}/fl.txt", fl_id=range(config["fl_params"]["N"])),
    output:
        "<output>/gp_map_04/navigability_per_fl_pop{pop_size}_vienna_RNA_gp_map.txt"
    resources:
        mem_mb_per_cpu=config["min_mem_per_cpu"],
        runtime=config["max_runtime"],
    wildcard_constraints:
        pop_size='[0-9]*',
        bp='[0-9]*'
    shell:
        "workflow/scripts/evo_simulations/compute_navigability_per_fl.py "
        "-w {input.walk_length_files} "  # assumes one per fitness landscape
        "-f {input.fitness_landscapes} "
        "-o {output} "
