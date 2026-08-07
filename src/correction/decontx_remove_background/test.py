import sys
import os
import pytest
import mudata as mu
from openpipeline_testutils.asserters import assert_annotation_objects_equal

## VIASH START
meta = {"resources_dir": "resources_test"}
## VIASH END

input_file = (
    f"{meta['resources_dir']}/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"
)

background_file = (
    f"{meta['resources_dir']}/pbmc_1k_protein_v3_raw_feature_bc_matrix.h5mu"
)


def assert_decontx_output(output_file):
    assert os.path.exists(output_file), "No output was created."

    input_mudata = mu.read_h5mu(input_file)
    output_mudata = mu.read_h5mu(output_file)

    assert_annotation_objects_equal(input_mudata.mod["prot"], output_mudata.mod["prot"])

    output_rna = output_mudata.mod["rna"]
    assert "decontx_counts" in output_rna.layers, (
        "Output should contain a 'decontx_counts' layer."
    )
    assert "decontx_contamination" in output_rna.obs, (
        "Output should contain a 'decontx_contamination' .obs column."
    )
    assert "decontx_clusters" in output_rna.obs, (
        "Output should contain a 'decontx_clusters' .obs column."
    )
    assert output_rna.obs["decontx_contamination"].between(0, 1).all(), (
        "Contamination values should be between 0 and 1."
    )


def test_execution_without_background(run_component, random_h5mu_path):
    output_file = random_h5mu_path()

    print("> Check whether decontx works without a background dataset")
    run_component(
        [
            "--input",
            input_file,
            "--output",
            output_file,
            "--output_compression",
            "gzip",
        ]
    )

    assert_decontx_output(output_file)


def test_execution_with_background(run_component, random_h5mu_path):
    output_file = random_h5mu_path()

    print("> Check whether decontx works with a background dataset")
    run_component(
        [
            "--input",
            input_file,
            "--background",
            background_file,
            "--output",
            output_file,
            "--output_compression",
            "gzip",
        ]
    )

    assert_decontx_output(output_file)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
