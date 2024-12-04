import mudata as mu
import sys
from pathlib import Path
import pytest

## VIASH START
meta = {
    "executable": "./target/executable/filter/filter_with_counts/filter_with_counts",
    "resources_dir": "resources_test/",
    "config": "/home/di/code/openpipeline/src/filter/filter_with_counts/config.vsh.yaml",
}

## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()


@pytest.fixture
def input_path():
    return f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"


@pytest.fixture
def input_h5mu(input_path):
    return mu.read_h5mu(input_path)


@pytest.fixture
def input_n_rna_obs(input_h5mu):
    return input_h5mu.mod["rna"].n_obs


@pytest.fixture
def input_n_prot_obs(input_h5mu):
    return input_h5mu.mod["prot"].n_obs


@pytest.fixture
def input_n_rna_vars(input_h5mu):
    return input_h5mu.mod["rna"].n_vars


@pytest.fixture
def input_n_prot_vars(input_h5mu):
    return input_h5mu.mod["prot"].n_vars


def test_filter_nothing(
    run_component,
    input_path,
    input_n_rna_obs,
    input_n_prot_obs,
    input_n_rna_vars,
    input_n_prot_vars,
):
    run_component(
        [
            "--input",
            input_path,
            "--output",
            "output-1.h5mu",
            "--min_cells_per_gene",
            "3",
            "--output_compression",
            "gzip",
        ]
    )
    assert Path("output-1.h5mu").is_file()
    mu_out = mu.read_h5mu("output-1.h5mu")
    assert "filter_with_counts" in mu_out.mod["rna"].obs
    assert "filter_with_counts" in mu_out.mod["rna"].var
    new_obs = mu_out.mod["rna"].n_obs
    new_vars = mu_out.mod["rna"].n_vars
    assert new_obs == input_n_rna_obs
    assert new_vars == input_n_rna_vars
    assert mu_out.mod["prot"].n_obs == input_n_prot_obs
    assert mu_out.mod["prot"].n_vars == input_n_prot_vars
    assert list(mu_out.mod["rna"].var["feature_types"].cat.categories) == [
        "Gene Expression"
    ]
    assert list(mu_out.mod["prot"].var["feature_types"].cat.categories) == [
        "Antibody Capture"
    ]


def test_filtering_a_little(
    run_component,
    input_path,
    input_n_rna_obs,
    input_n_prot_obs,
    input_n_rna_vars,
    input_n_prot_vars,
):
    run_component(
        [
            "--input",
            input_path,
            "--output",
            "output-2.h5mu",
            "--modality",
            "rna",
            "--min_counts",
            "200",
            "--max_counts",
            "5000000",
            "--min_genes_per_cell",
            "200",
            "--max_genes_per_cell",
            "1500000",
            "--min_cells_per_gene",
            "10",
            "--do_subset",
        ]
    )
    assert Path("output-2.h5mu").is_file()
    mu_out = mu.read_h5mu("output-2.h5mu")
    new_obs = mu_out.mod["rna"].n_obs
    new_vars = mu_out.mod["rna"].n_vars
    assert new_obs < input_n_rna_obs
    assert new_vars < input_n_rna_vars
    assert mu_out.mod["prot"].n_obs == input_n_prot_obs
    assert mu_out.mod["prot"].n_vars == input_n_prot_vars
    assert list(mu_out.mod["rna"].var["feature_types"].cat.categories) == [
        "Gene Expression"
    ]
    assert list(mu_out.mod["prot"].var["feature_types"].cat.categories) == [
        "Antibody Capture"
    ]


def test_filter_cells_without_counts(run_component, input_h5mu, tmp_path):
    # create_an_empty_cell
    obs_to_remove = input_h5mu.mod["rna"].obs.index[0]
    input_h5mu.mod["rna"].X[0] = 0
    temp_h5mu_path = tmp_path / "temp.h5mu"
    input_h5mu.write(temp_h5mu_path)
    run_component(
        [
            "--input",
            temp_h5mu_path,
            "--output",
            "output-3.h5mu",
            "--min_cells_per_gene",
            "0",
        ]
    )
    assert Path("output-3.h5mu").is_file()
    mu_out = mu.read_h5mu("output-3.h5mu")
    assert mu_out.mod["rna"].obs.at[obs_to_remove, "filter_with_counts"] is False
    assert "mitochondrial" not in mu_out.mod["rna"].var


def test_filter_using_different_layer(
    run_component,
    input_h5mu,
    tmp_path,
    input_n_rna_obs,
    input_n_prot_obs,
    input_n_rna_vars,
    input_n_prot_vars,
):
    # move X to different input layer
    input_h5mu.mod["rna"].layers["test_layer"] = input_h5mu.mod["rna"].X.copy()
    input_h5mu.mod["rna"].X = None

    temp_h5mu_path = tmp_path / "temp.h5mu"
    input_h5mu.write(temp_h5mu_path)
    run_component(
        [
            "--input",
            temp_h5mu_path,
            "--output",
            "output-4.h5mu",
            "--modality",
            "rna",
            "--min_counts",
            "200",
            "--max_counts",
            "5000000",
            "--min_genes_per_cell",
            "200",
            "--max_genes_per_cell",
            "1500000",
            "--min_cells_per_gene",
            "10",
            "--layer",
            "test_layer",
            "--do_subset",
        ]
    )
    assert Path("output-4.h5mu").is_file()
    mu_out = mu.read_h5mu("output-2.h5mu")
    new_obs = mu_out.mod["rna"].n_obs
    new_vars = mu_out.mod["rna"].n_vars
    assert new_obs < input_n_rna_obs
    assert new_vars < input_n_rna_vars
    assert mu_out.mod["prot"].n_obs == input_n_prot_obs
    assert mu_out.mod["prot"].n_vars == input_n_prot_vars
    assert list(mu_out.mod["rna"].var["feature_types"].cat.categories) == [
        "Gene Expression"
    ]
    assert list(mu_out.mod["prot"].var["feature_types"].cat.categories) == [
        "Antibody Capture"
    ]


if __name__ == "__main__":
    exit(pytest.main([__file__]))
