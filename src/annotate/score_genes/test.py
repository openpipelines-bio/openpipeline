import mudata as mu

input_file = f"{meta['resources_dir']}/merge_test_data/pbmc_1k_protein_v3_filtered_feature_bc_matrix_rna.h5mu"
output_file = f"{meta['temp_dir']}/output.h5mu"

# generate gene list file
gene_list = ["MCM5", "PCNA", "TYMS"]
gene_list_file = f"{meta['temp_dir']}/gene_list.txt"

with open(gene_list_file, 'w') as file:
    for gene in gene_list:
        file.write(gene + '\n')


def test_cell_scoring(run_component):

    run_component([
        "--input", input_file,
        "--gene_list", gene_list_file,
        "--output", output_file,
        "--score_name", 'cell_cycle_score'
    ])

    output = mu.read(output_file)
    
    # check output
    expected_rna_obs_cols = ["cell_cycle_score"]
    for col in expected_rna_obs_cols:
        assert col in output.mod["rna"].obs.columns, f"could not find columns .mod['rna'].obs['{col}']"