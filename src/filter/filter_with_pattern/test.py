import mudata as mu
import sys
from pathlib import Path
import pytest
import numpy as np

## VIASH START
meta = {
    "executable": "./target/executable/filter/filter_with_pattern/filter_with_pattern",
    "resources_dir": "resources_test/",
    "config": "/home/di/code/openpipeline/src/filter/filter_with_pattern/config.vsh.yaml",
}

## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()


@pytest.fixture
def input_path():
    return f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"


def test_simple_execution(run_component, input_path, tmp_path):
    output_path = tmp_path / "output.h5mu"
    run_component(
        [
            "--input",
            input_path,
            "--output",
            output_path,
            "--var_gene_names",
            "gene_symbol",
            "--pattern",
            "MIR\\d+",
            "--pattern",
            "MIR\\d+",
            "--pattern",
            "LINC\\d+",
            "--pattern",
            "AC\\d+",
            "--pattern",
            "AP\\d+",
            "--do_subset",
            "True",
            "--output_compression",
            "gzip",
        ]
    )
    assert Path(output_path).is_file()
    mu_out = mu.read_h5ad(output_path, mod="rna")
    mu_in = mu.read_h5ad(input_path, mod="rna")

    var_out = mu_out.var

    assert "filter_with_pattern" in var_out, (
        "Expected column with filter results to be present"
    )
    assert var_out["filter_with_pattern"].dtype == np.dtype("bool"), (
        "Expected filter column to be boolean"
    )

    mu_in.shape[0] == mu_out.shape[0], "Number of observations should be unaltered"
    mu_out.shape[1] < mu_in.shape[1], "Genes should have been filtered"


def test_filter_without_subset(run_component, input_path, tmp_path):
    output_path = tmp_path / "output.h5mu"
    run_component(
        [
            "--input",
            input_path,
            "--output",
            output_path,
            "--var_gene_names",
            "gene_symbol",
            "--pattern",
            "MIR\\d+",
            "--pattern",
            "MIR\\d+",
            "--pattern",
            "LINC\\d+",
            "--pattern",
            "AC\\d+",
            "--pattern",
            "AP\\d+",
            "--do_subset",
            "False",
            "--output_compression",
            "gzip",
        ]
    )
    assert Path(output_path).is_file()
    mu_out = mu.read_h5ad(output_path, mod="rna")
    mu_in = mu.read_h5ad(input_path, mod="rna")

    var_out = mu_out.var

    assert "filter_with_pattern" in var_out, (
        "Expected column with filter results to be present"
    )
    assert var_out["filter_with_pattern"].dtype == np.dtype("bool"), (
        "Expected filter column to be boolean"
    )

    mu_in.shape[0] == mu_out.shape[0], "Number of observations should be unaltered"
    mu_out.shape[1] == mu_in.shape[1], "Genes should not have been filtered"


if __name__ == "__main__":
    exit(pytest.main([__file__]))
