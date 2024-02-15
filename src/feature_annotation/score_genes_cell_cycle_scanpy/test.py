import pytest
import sys
import mudata as mu

input_file = f"{meta['resources_dir']}/merge_test_data/pbmc_1k_protein_v3_filtered_feature_bc_matrix_rna.h5mu"


@pytest.fixture
def s_genes(tmp_path):
    result = tmp_path / "s_genes.txt"
    s_genes = ["MCM5", "PCNA", "TYMS"]
    with result.open('w') as open_s_genes_file:
        open_s_genes_file.write("\n".join(s_genes))
    return result


@pytest.fixture
def g2m_genes(tmp_path):
    result = tmp_path / "s_genes.txt"
    g2m_genes = ["UBE2C", "BIRC5", "TPX2"]
    with result.open('w') as open_g2m_genes_file:
        open_g2m_genes_file.write("\n".join(g2m_genes))
    return result


def test_cell_scoring_cell_cycle(run_component, tmp_path, s_genes, g2m_genes):

    output_file = tmp_path / "output.h5mu"

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


if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))