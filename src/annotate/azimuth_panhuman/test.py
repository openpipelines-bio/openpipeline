import os
import re
import subprocess
import sys

import mudata as mu
import numpy as np
import pandas as pd
import pytest
from openpipeline_testutils.asserters import assert_annotation_objects_equal

## VIASH START
meta = {
    "executable": "./target/docker/annotate/azimuth_panhuman/azimuth_panhuman",
    "resources_dir": "resources_test/",
    "cpus": 4,
    "memory_gb": 20,
    "config": "src/annotate/azimuth_panhuman/config.vsh.yaml",
}
## VIASH END

# Raw (unnormalized) counts, exactly as produced by cellranger -> h5mu
# conversion: Azimuth performs its own normalization and expects integer
# counts, so this is used as-is (no log-normalization fixture needed).
input_file = (
    f"{meta['resources_dir']}/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"
)
# Pre-downloaded v1 model, identical to the one panhumanpy downloads itself.
model_file = f"{meta['resources_dir']}/panhumanpy_inference_model_v1.keras"


def test_simple_execution(run_component, random_h5mu_path):
    output_file = random_h5mu_path()

    run_component(
        [
            "--input",
            input_file,
            "--input_var_gene_names",
            "gene_symbol",
            "--output",
            output_file,
        ]
    )

    assert os.path.exists(output_file), "Output file does not exist"

    input_mudata = mu.read_h5mu(input_file)
    output_mudata = mu.read_h5mu(output_file)

    assert_annotation_objects_equal(input_mudata.mod["prot"], output_mudata.mod["prot"])

    output_rna = output_mudata.mod["rna"]

    # Only the parametrized columns are added to .obs (the input has none):
    # the final-level prediction/confidence and the refined broad/medium/fine
    # labels, under their --output_obs_* default names.
    expected_obs_cols = {
        "azimuth_pred",
        "azimuth_probability",
        "azimuth_broad",
        "azimuth_medium",
        "azimuth_fine",
    }
    assert set(output_rna.obs.columns) == expected_obs_cols

    predictions = output_rna.obs["azimuth_broad"]
    assert not all(predictions.isna()), "Not all predictions should be NA"

    confidence = output_rna.obs["azimuth_probability"]
    assert all(0 <= value <= 1 for value in confidence), (
        ".obs at azimuth_probability has values outside the range [0, 1]"
    )

    assert "X_azimuth" in output_rna.obsm, "Embeddings were not stored in .obsm"
    assert "X_azimuth_umap" in output_rna.obsm, "UMAP was not stored in .obsm"
    assert output_rna.obsm["X_azimuth"].shape[0] == output_rna.n_obs
    assert output_rna.obsm["X_azimuth_umap"].shape == (output_rna.n_obs, 2)


def test_no_refine_labels_and_no_embeddings(run_component, random_h5mu_path):
    output_file = random_h5mu_path()

    run_component(
        [
            "--input",
            input_file,
            "--input_var_gene_names",
            "gene_symbol",
            "--refine_labels",
            "false",
            "--extract_embeddings",
            "false",
            "--umap_embeddings",
            "false",
            "--output",
            output_file,
        ]
    )

    assert os.path.exists(output_file), "Output file does not exist"

    output_rna = mu.read_h5mu(output_file).mod["rna"]

    # refine_labels=false: no broad/medium/fine columns
    assert set(output_rna.obs.columns) == {"azimuth_pred", "azimuth_probability"}

    # extract_embeddings=false / umap_embeddings=false: no obsm entries added
    assert "X_azimuth" not in output_rna.obsm
    assert "X_azimuth_umap" not in output_rna.obsm


def test_duplicate_categorical_index_entry(
    run_component, random_h5mu_path, write_mudata_to_file
):
    """Component should not raise on duplicate .var index entries."""
    output_file = random_h5mu_path()

    input_mdata = mu.read_h5mu(input_file)
    indices = input_mdata.mod["rna"].var["gene_symbol"].to_list()
    indices[0] = indices[1]  # duplicate the first index entry
    input_mdata.mod["rna"].var["dup_idx"] = pd.CategoricalIndex(indices)

    dup_input_file = write_mudata_to_file(input_mdata)

    run_component(
        [
            "--input",
            dup_input_file,
            "--input_var_gene_names",
            "dup_idx",
            "--output",
            output_file,
        ]
    )

    assert os.path.exists(output_file), "Output file does not exist"


def test_model_version_v0_and_custom_outputs(run_component, random_h5mu_path):
    output_file = random_h5mu_path()

    run_component(
        [
            "--input",
            input_file,
            "--input_var_gene_names",
            "gene_symbol",
            "--model_version",
            "v0",
            "--umap_n_neighbors",
            "10",
            "--umap_metric",
            "euclidean",
            "--output_obs_predictions",
            "my_pred",
            "--output_obs_probability",
            "my_probability",
            "--output_obs_predictions_broad",
            "my_broad",
            "--output_obs_predictions_medium",
            "my_medium",
            "--output_obs_predictions_fine",
            "my_fine",
            "--output_obsm_embedding",
            "my_embedding",
            "--output_obsm_umap",
            "my_umap",
            "--output_compression",
            "gzip",
            "--output",
            output_file,
        ]
    )

    assert os.path.exists(output_file), "Output file does not exist"

    output_rna = mu.read_h5mu(output_file).mod["rna"]
    assert "my_pred" in output_rna.obs, (
        "Predictions were not stored under the custom obs key"
    )
    assert "my_probability" in output_rna.obs, (
        "Probability was not stored under the custom obs key"
    )
    for custom in ["my_broad", "my_medium", "my_fine"]:
        assert custom in output_rna.obs, f"{custom} missing from .obs"
    # No columns stored under the default names
    assert set(output_rna.obs.columns) == {
        "my_pred",
        "my_probability",
        "my_broad",
        "my_medium",
        "my_fine",
    }

    assert "my_embedding" in output_rna.obsm, (
        "Embeddings were not stored under the custom obsm key"
    )
    assert "my_umap" in output_rna.obsm, "UMAP was not stored under the custom obsm key"
    assert "X_azimuth" not in output_rna.obsm
    assert "X_azimuth_umap" not in output_rna.obsm


def test_provided_model(run_component, random_h5mu_path):
    output_file = random_h5mu_path()

    component_output = run_component(
        [
            "--input",
            input_file,
            "--input_var_gene_names",
            "gene_symbol",
            "--model",
            model_file,
            "--model_version",
            "v1",
            "--output",
            output_file,
        ]
    ).decode("utf-8")

    # The provided model is loaded instead of being downloaded by panhumanpy
    assert "Using provided model" in component_output
    assert "Downloading model" not in component_output

    assert os.path.exists(output_file), "Output file does not exist"

    input_mudata = mu.read_h5mu(input_file)
    output_mudata = mu.read_h5mu(output_file)

    assert_annotation_objects_equal(input_mudata.mod["prot"], output_mudata.mod["prot"])

    output_rna = output_mudata.mod["rna"]

    assert set(output_rna.obs.columns) == {
        "azimuth_pred",
        "azimuth_probability",
        "azimuth_broad",
        "azimuth_medium",
        "azimuth_fine",
    }

    predictions = output_rna.obs["azimuth_broad"]
    assert not all(predictions.isna()), "Not all predictions should be NA"

    confidence = output_rna.obs["azimuth_probability"]
    assert all(0 <= value <= 1 for value in confidence), (
        ".obs at azimuth_probability has values outside the range [0, 1]"
    )

    assert output_rna.obsm["X_azimuth"].shape[0] == output_rna.n_obs
    assert output_rna.obsm["X_azimuth_umap"].shape == (output_rna.n_obs, 2)


def test_fail_normalized_input(run_component, random_h5mu_path, write_mudata_to_file):
    output_file = random_h5mu_path()

    input_mdata = mu.read_h5mu(input_file)
    adata = input_mdata.mod["rna"]
    X = adata.X.astype(float).tocsr()
    row_sums = np.asarray(X.sum(axis=1)).ravel()
    row_sums[row_sums == 0] = 1
    X = X.multiply(10000 / row_sums[:, None]).tocsr()
    X.data = np.log1p(X.data)
    adata.X = X

    lognorm_input_file = write_mudata_to_file(input_mdata)

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                lognorm_input_file,
                "--input_var_gene_names",
                "gene_symbol",
                "--output",
                output_file,
            ]
        )
    assert re.search(
        r"detected non-integer values.*already normalized",
        err.value.stdout.decode("utf-8"),
    )
    assert not os.path.exists(output_file), "Output file should not have been created"


def test_fail_insufficient_gene_overlap(run_component, random_h5mu_path):
    output_file = random_h5mu_path()

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_file,
                "--input_var_gene_names",
                "gene_symbol",
                "--input_reference_gene_overlap",
                "1000000",
                "--output",
                output_file,
            ]
        )
    assert re.search(
        r"intersection of genes between the query and reference dataset is too small",
        err.value.stdout.decode("utf-8"),
    )
    assert not os.path.exists(output_file), "Output file should not have been created"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
