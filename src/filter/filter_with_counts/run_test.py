import mudata as mu
import logging
from sys import stdout, exit
from pathlib import Path
import pytest

## VIASH START
meta = {
    'executable': './target/docker/filter/filter_with_counts/filter_with_counts',
    'resources_dir': 'resources_test/',
    'config': "/home/di/code/openpipeline/src/filter/filter_with_counts/config.vsh.yaml"
}

## VIASH END

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)


@pytest.fixture
def input_path():
    return f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"

@pytest.fixture
def input_h5mu(input_path):
    return mu.read_h5mu(input_path)

@pytest.fixture
def input_n_rna_obs(input_h5mu):
    return input_h5mu.mod['rna'].n_obs

@pytest.fixture
def input_n_prot_obs(input_h5mu):
    return input_h5mu.mod['prot'].n_obs

@pytest.fixture
def input_n_rna_vars(input_h5mu):
    return input_h5mu.mod['rna'].n_vars

@pytest.fixture
def input_n_prot_vars(input_h5mu):
    return input_h5mu.mod['prot'].n_vars

def test_filter_nothing(run_component, input_path,
                        input_n_rna_obs, input_n_prot_obs,
                        input_n_rna_vars, input_n_prot_vars):
    run_component([
        "--input", input_path,
        "--output", "output-1.h5mu",
        "--min_cells_per_gene", "3"
        ])
    assert Path("output-1.h5mu").is_file()
    mu_out = mu.read_h5mu("output-1.h5mu")
    assert "filter_with_counts" in mu_out.mod["rna"].obs
    assert "filter_with_counts" in mu_out.mod["rna"].var
    new_obs = mu_out.mod['rna'].n_obs
    new_vars = mu_out.mod['rna'].n_vars
    assert new_obs == input_n_rna_obs
    assert new_vars == input_n_rna_vars
    assert mu_out.mod['prot'].n_obs == input_n_prot_obs
    assert mu_out.mod['prot'].n_vars == input_n_prot_vars
    assert list(mu_out.mod['rna'].var['feature_types'].cat.categories) == ["Gene Expression"]
    assert list(mu_out.mod['prot'].var['feature_types'].cat.categories) == ["Antibody Capture"]

def test_filtering_a_little(run_component, input_path,
                            input_n_rna_obs, input_n_prot_obs,
                            input_n_rna_vars, input_n_prot_vars):
    run_component([
        "--input", input_path,
        "--output", "output-2.h5mu",
        "--modality", "rna",
        "--min_counts", "200",
        "--max_counts", "5000000",
        "--min_genes_per_cell", "200",
        "--max_genes_per_cell", "1500000",
        "--min_cells_per_gene", "10",
        "--min_fraction_mito", "0",
        "--max_fraction_mito", "0.2",
        "--var_gene_names", "gene_symbol",
        "--do_subset"])
    assert Path("output-2.h5mu").is_file()
    mu_out = mu.read_h5mu("output-2.h5mu")
    new_obs = mu_out.mod['rna'].n_obs
    new_vars = mu_out.mod['rna'].n_vars
    assert new_obs < input_n_rna_obs
    assert new_vars < input_n_rna_vars
    assert mu_out.mod['prot'].n_obs == input_n_prot_obs
    assert mu_out.mod['prot'].n_vars == input_n_prot_vars
    assert list(mu_out.mod['rna'].var['feature_types'].cat.categories) == ["Gene Expression"]
    assert list(mu_out.mod['prot'].var['feature_types'].cat.categories) == ["Antibody Capture"]

def test_filter_cells_without_counts(run_component, input_h5mu, tmp_path):
    # create_an_empty_cell
    obs_to_remove = input_h5mu.mod['rna'].obs.index[0]
    input_h5mu.mod['rna'].X[0] = 0
    temp_h5mu_path = tmp_path / "temp.h5mu"
    input_h5mu.write(temp_h5mu_path)
    run_component([
        "--input", temp_h5mu_path,
        "--output", "output-3.h5mu",
        "--min_cells_per_gene", "0",
    ])
    assert Path("output-3.h5mu").is_file()
    mu_out = mu.read_h5mu("output-3.h5mu")
    assert mu_out.mod['rna'].obs.at[obs_to_remove, 'filter_with_counts'] == False

def test_filter_mitochondrial(run_component, input_path,
                              input_n_rna_obs, input_n_prot_obs,
                              input_n_rna_vars, input_n_prot_vars):
    run_component([
        "--input", input_path,
        "--output", "output-4.h5mu",
        "--var_gene_names", "gene_symbol",
        "--max_fraction_mito", "0.2",
        "--do_subset"
        ])
    assert Path("output-4.h5mu").is_file()
    mu_out = mu.read_h5mu("output-4.h5mu")
    new_obs = mu_out.mod['rna'].n_obs
    new_vars = mu_out.mod['rna'].n_vars
    assert new_obs < input_n_rna_obs
    assert new_vars == input_n_rna_vars
    assert mu_out.mod['prot'].n_obs == input_n_prot_obs
    assert mu_out.mod['prot'].n_vars == input_n_prot_vars
    assert list(mu_out.mod['rna'].var['feature_types'].cat.categories) == ["Gene Expression"]
    assert list(mu_out.mod['prot'].var['feature_types'].cat.categories) == ["Antibody Capture"]

def test_filter_mitochondrial_regex(run_component, input_path,
                                    input_n_rna_obs, input_n_prot_obs,
                                    input_n_rna_vars, input_n_prot_vars):
    run_component([
        "--input", input_path,
        "--output", "output-5.h5mu",
        "--var_gene_names", "gene_symbol",
        "--max_fraction_mito", "0.2",
        "--mitochondrial_gene_regex", "^[M][T]-",
        "--do_subset"
        ])
    assert Path("output-5.h5mu").is_file()
    mu_out = mu.read_h5mu("output-5.h5mu")
    new_obs = mu_out.mod['rna'].n_obs
    new_vars = mu_out.mod['rna'].n_vars
    assert new_obs < input_n_rna_obs
    assert new_vars == input_n_rna_vars
    assert mu_out.mod['prot'].n_obs == input_n_prot_obs
    assert mu_out.mod['prot'].n_vars == input_n_prot_vars
    assert list(mu_out.mod['rna'].var['feature_types'].cat.categories) == ["Gene Expression"]
    assert list(mu_out.mod['prot'].var['feature_types'].cat.categories) == ["Antibody Capture"]

def test_filter_mitochondrial_column_not_set(run_component, input_path,
                                             input_n_rna_obs, input_n_prot_obs,
                                             input_n_rna_vars, input_n_prot_vars):
    run_component([
        "--input", input_path,
        "--output", "output-6.h5mu",
        "--max_fraction_mito", "0.2",
        "--mitochondrial_gene_regex", "^[mM][tT]-",
        "--do_subset"
        ])
    assert Path("output-6.h5mu").is_file()
    mu_out = mu.read_h5mu("output-6.h5mu")
    new_obs = mu_out.mod['rna'].n_obs
    new_vars = mu_out.mod['rna'].n_vars
    assert new_obs == input_n_rna_obs # Stays the same because filtering on the wrong column
    assert new_vars == input_n_rna_vars
    assert mu_out.mod['prot'].n_obs == input_n_prot_obs
    assert mu_out.mod['prot'].n_vars == input_n_prot_vars
    assert list(mu_out.mod['rna'].var['feature_types'].cat.categories) == ["Gene Expression"]
    assert list(mu_out.mod['prot'].var['feature_types'].cat.categories) == ["Antibody Capture"]


if __name__ == "__main__":
    exit(pytest.main([__file__], plugins=["viashpy"]))
