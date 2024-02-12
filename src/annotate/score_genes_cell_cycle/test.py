import mudata as mu

## VIASH START
meta = {
    "resources_dir": "resources_test"
}
## VIASH END

input_file = f"{meta['resources_dir']}/merge_test_data/pbmc_1k_protein_v3_filtered_feature_bc_matrix_rna.h5mu"
output_file = f"{meta['temp_dir']}/output.h5mu"

# generate gene list files
s_genes = ["MCM5", "PCNA", "TYMS"]
s_genes_file = f"{meta['temp_dir']}/s_genes.txt"

g2m_genes = ["UBE2C", "BIRC5", "TPX2"]
g2m_genes_file = f"{meta['temp_dir']}/g2m_genes.txt"


with open(s_genes_file, 'w') as file:
    for gene in s_genes:
        file.write(gene + '\n')

with open(g2m_genes_file, 'w') as file:
    for gene in g2m_genes:
        file.write(gene + '\n')


def test_cell_scoring_cell_cycle(run_component):

    run_component([
        "--input", input_file,
        "--s_genes", s_genes,
        "--g2m_genes", g2m_genes,
        "--output", output_file
    ])

    output = mu.read(output_file)

    # check output
    expected_rna_obs_cols = ["S_score", "G2M_score", 'phase']
    for col in expected_rna_obs_cols:
        assert col in output.mod["rna"].obs.columns, f"could not find columns .mod['rna'].obs['{col}']"