import os
import re
import subprocess
import sys

import mudata as mu
import pytest
from openpipeline_testutils.asserters import assert_annotation_objects_equal

## VIASH START
meta = {
    "executable": "./target/docker/annotate/azimuth_seurat/azimuth_seurat",
    "resources_dir": "resources_test/",
    "cpus": 4,
    "memory_gb": 20,
    "config": "src/annotate/azimuth_seurat/config.vsh.yaml",
}
## VIASH END

# Raw (unnormalized) counts: RunAzimuth() recomputes SCT residuals from raw
# counts internally, so this is used as-is (no pre-normalization needed).
input_file = (
    f"{meta['resources_dir']}/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"
)


def test_simple_execution(run_component, random_h5mu_path):
    output_file = random_h5mu_path()

    run_component(
        [
            "--input",
            input_file,
            "--reference",
            "pbmcref",
            "--output",
            output_file,
        ]
    )

    assert os.path.exists(output_file), "Output file does not exist"

    input_mudata = mu.read_h5mu(input_file)
    output_mudata = mu.read_h5mu(output_file)

    assert_annotation_objects_equal(input_mudata.mod["prot"], output_mudata.mod["prot"])

    output_rna = output_mudata.mod["rna"]

    predicted_cols = [c for c in output_rna.obs.columns if c.startswith("predicted.")]
    assert predicted_cols, "No 'predicted.*' annotation columns were added"

    assert "mapping.score" in output_rna.obs
    assert all(0 <= value <= 1 for value in output_rna.obs["mapping.score"]), (
        ".obs at mapping.score has values outside the range [0, 1]"
    )

    assert "ref.umap" in output_rna.obsm, "UMAP was not projected onto the reference"
    assert output_rna.obsm["ref.umap"].shape == (output_rna.n_obs, 2)


def test_annotation_levels_subset(run_component, random_h5mu_path):
    output_file = random_h5mu_path()

    run_component(
        [
            "--input",
            input_file,
            "--reference",
            "pbmcref",
            "--annotation_levels",
            "celltype.l1",
            "--output",
            output_file,
        ]
    )

    assert os.path.exists(output_file), "Output file does not exist"

    output_rna = mu.read_h5mu(output_file).mod["rna"]
    assert "predicted.celltype.l1" in output_rna.obs
    assert "predicted.celltype.l2" not in output_rna.obs


def test_fail_invalid_reference(run_component, random_h5mu_path):
    output_file = random_h5mu_path()

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_file,
                "--reference",
                "not_a_real_reference",
                "--output",
                output_file,
            ]
        )
    assert re.search(
        r"Could not find a reference for",
        err.value.stdout.decode("utf-8"),
    )
    assert not os.path.exists(output_file), "Output file should not have been created"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
