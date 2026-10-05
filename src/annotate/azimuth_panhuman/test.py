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

    # Note: "level_zero_labels" is replaced by "azimuth_broad" once label
    # refinement runs (the default), rather than existing alongside it.
    # "final_level_labels"/"final_level_confidence" are renamed to the
    # --output_obs_predictions/--output_obs_probability defaults.
    expected_obs_cols = {
        "full_hierarchical_labels",
        "azimuth_pred",
        "azimuth_probability",
        "azimuth_broad",
        "azimuth_medium",
        "azimuth_fine",
    }
    assert expected_obs_cols.issubset(output_rna.obs.keys()), (
        f"Missing expected .obs columns: {expected_obs_cols - set(output_rna.obs.keys())}"
    )

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


def test_no_refine_detailed_output_and_cell_ontology(run_component, random_h5mu_path):
    output_file = random_h5mu_path()

    run_component(
        [
            "--input",
            input_file,
            "--input_var_gene_names",
            "gene_symbol",
            "--refine_labels",
            "false",
            "--output_mode",
            "detailed",
            "--map_to_cl",
            "level_zero_labels",
            "--include_cl_id",
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

    # refine_labels=false: no azimuth_broad/medium/fine columns
    assert "azimuth_broad" not in output_rna.obs
    assert "azimuth_medium" not in output_rna.obs
    assert "azimuth_fine" not in output_rna.obs

    # output_mode=detailed: per-level label columns are present
    assert "level_1_labels" in output_rna.obs

    # map_to_cl + include_cl_id: CL columns derived from level_zero_labels
    assert "level_zero_labels_CL" in output_rna.obs
    assert "level_zero_labels_CL_ID" in output_rna.obs

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
    for custom, original in [
        ("my_broad", "azimuth_broad"),
        ("my_medium", "azimuth_medium"),
        ("my_fine", "azimuth_fine"),
    ]:
        assert custom in output_rna.obs, f"{custom} missing from .obs"
        assert original not in output_rna.obs, f"{original} should be renamed"
    assert "final_level_labels" not in output_rna.obs
    assert "final_level_confidence" not in output_rna.obs

    assert "my_embedding" in output_rna.obsm, (
        "Embeddings were not stored under the custom obsm key"
    )
    assert "my_umap" in output_rna.obsm, "UMAP was not stored under the custom obsm key"
    assert "X_azimuth" not in output_rna.obsm
    assert "X_azimuth_umap" not in output_rna.obsm


def test_provided_model_is_used_instead_of_download(
    run_component, random_h5mu_path, tmp_path
):
    # A corrupt model can only make the run fail if the provided file is
    # actually loaded; had panhumanpy downloaded the weights instead, the
    # run would have succeeded.
    fake_model = tmp_path / "fake_model.keras"
    fake_model.write_text("not a keras model")

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_file,
                "--input_var_gene_names",
                "gene_symbol",
                "--model",
                str(fake_model),
                "--output",
                random_h5mu_path(),
            ]
        )
    assert "Using provided model" in err.value.stdout.decode("utf-8")


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
