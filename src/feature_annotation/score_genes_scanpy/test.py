import pytest
import sys
import mudata as mu

input_file = f"{meta['resources_dir']}/merge_test_data/pbmc_1k_protein_v3_filtered_feature_bc_matrix_rna.h5mu"


@pytest.fixture
def gene_list_file(tmp_path):
    result = tmp_path / "s_genes.txt"
    gene_list = ["MCM5", "PCNA", "TYMS"]
    with result.open('w') as open_gene_list_file:
        open_gene_list_file.write("\n".join(gene_list))
    return result


def test_cell_scoring(run_component, tmp_path, gene_list_file):
    
    output_file = tmp_path / "output.h5mu"

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


if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))