import sys
import tskit


if __name__ == '__main__':
    program_name = sys.argv[0]
    program_arguments = sys.argv[1:]
    population_id = int(program_arguments[0])
    input_ts_file = program_arguments[1]
    output_al_file = program_arguments[2]
    number_of_reference_populations = int(program_arguments[3])
    offset_number = int(program_arguments[4])
    reference_population_ids = {}
    for pid in range(number_of_reference_populations):
        reference_population_ids[program_arguments[5 + pid]] = pid
    ts = tskit.load(input_ts_file)
    ts = ts.simplify(keep_input_roots=True)
    sample_node_ids = ts.samples(population=population_id, time=0)
    with open(output_al_file, "w") as file_output:
        file_output.write(f"#Population_label_and_ID:")
        for plbl, pid in reference_population_ids.items():
            file_output.write(f"\t{plbl}={pid}")
        file_output.write(f"\n#Individual_Haplotype_Node_ID\tPhysical_Position_Start(Inclusive)\tPhysical_Position_End(Exclusive)\tAncestral_Label_Index\n")
    is_first_haplotype = True
    for query_node_id_index, query_node_id in enumerate(sample_node_ids):
        st_prev = -1
        ed_prev = -1
        av_prev = -1
        with open(output_al_file, "a") as file_output:
            if is_first_haplotype:
                file_output.write(f"ADMIX_{query_node_id_index // 2}_0")
                is_first_haplotype = False
            else:
                file_output.write(f"ADMIX_{query_node_id_index // 2}_1")
                is_first_haplotype = True
            for tree in ts.trees():
                st = int(tree.interval[0])
                ed = int(tree.interval[1])
                node_id = query_node_id
                while tree.parent(node_id) != tskit.NULL:
                    node_id = tree.parent(node_id)
                av = tree.population(node_id)
                av -= offset_number
                if av_prev == av:
                    ed_prev = ed
                else:
                    if av_prev != -1:
                        file_output.write(f"\t{st_prev}\t{ed_prev}\t{av_prev}")
                    st_prev = st
                    ed_prev = ed
                    av_prev = av
            file_output.write(f"\t{st_prev}\t{ed_prev}\t{av_prev}\n")
    sys.exit(0)
