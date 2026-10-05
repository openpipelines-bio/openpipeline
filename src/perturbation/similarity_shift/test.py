import re
import sys
from subprocess import CalledProcessError

import anndata as ad
import mudata as mu
import numpy as np
import pandas as pd
import pytest

## VIASH START
meta = {
    "executable": "./target/executable/perturbation/similarity_shift/similarity_shift",
    "resources_dir": "src/utils",
    "config": "./src/perturbation/similarity_shift/config.vsh.yaml",
}
## VIASH END

HEALTHY = np.array([1.0, 0.0, 0.0])
DISEASE = np.array([0.0, 1.0, 0.0])


@pytest.fixture
def original_path(random_h5mu_path):
    rna = ad.AnnData(
        X=np.zeros((2, 1)),
        obs=pd.DataFrame(index=["cell_a", "cell_b"]),
    )
    # cell_a sits on the disease centroid, cell_b halfway
    rna.obsm["X_geneformer"] = np.array([[0.0, 2.0, 0.0], [1.0, 1.0, 0.0]])
    rna.uns["perturbation_centroids"] = pd.DataFrame(
        [HEALTHY, DISEASE], index=["healthy", "disease"], columns=["0", "1", "2"]
    )
    path = random_h5mu_path()
    mu.MuData({"rna": rna}).write_h5mu(path)
    return path


@pytest.fixture
def virtual_path(random_h5mu_path):
    obs = pd.DataFrame(
        {
            "perturbation_cell_id": ["cell_a", "cell_a", "cell_b"],
            "perturbation_gene_id": ["ENSG1", "ENSG2", "ENSG1"],
            "perturbation_gene_name": ["G1", "G2", "G1"],
        },
        index=["cell_a@ENSG1", "cell_a@ENSG2", "cell_b@ENSG1"],
    )
    rna = ad.AnnData(X=np.zeros((3, 0)), obs=obs)
    rna.obsm["X_geneformer"] = np.array(
        [
            [3.0, 0.0, 0.0],  # cell_a moved onto the healthy centroid
            [0.0, 5.0, 0.0],  # cell_a unchanged in direction
            [0.0, 0.0, 1.0],  # cell_b moved away from both
        ]
    )
    path = random_h5mu_path()
    mu.MuData({"rna": rna}).write_h5mu(path)
    return path


def test_shift_values(run_component, virtual_path, original_path, tmp_path):
    output = tmp_path / "shift.csv"
    run_component(
        ["--input", virtual_path, "--original", original_path, "--output", output]
    )
    result = pd.read_csv(output)
    assert list(result.columns) == [
        "cell_id",
        "gene_id",
        "gene_name",
        "cos_original_healthy",
        "cos_original_disease",
        "cos_perturbed_healthy",
        "cos_perturbed_disease",
        "shift_healthy",
        "shift_disease",
    ]
    assert list(result["cell_id"]) == ["cell_a", "cell_a", "cell_b"]
    assert list(result["gene_name"]) == ["G1", "G2", "G1"]

    half = 1 / np.sqrt(2)
    assert np.allclose(result["shift_healthy"], [1.0, 0.0, -half])
    assert np.allclose(result["shift_disease"], [-1.0, 0.0, -half])
    assert np.allclose(result["cos_original_healthy"], [0.0, 0.0, half])


def test_missing_original_cell_errors(
    run_component, virtual_path, random_h5mu_path, tmp_path
):
    rna = ad.AnnData(X=np.zeros((1, 1)), obs=pd.DataFrame(index=["cell_a"]))
    rna.obsm["X_geneformer"] = np.ones((1, 3))
    rna.uns["perturbation_centroids"] = pd.DataFrame(
        [HEALTHY, DISEASE], index=["healthy", "disease"], columns=["0", "1", "2"]
    )
    path = random_h5mu_path()
    mu.MuData({"rna": rna}).write_h5mu(path)
    with pytest.raises(CalledProcessError) as err:
        run_component(
            [
                "--input",
                virtual_path,
                "--original",
                path,
                "--output",
                tmp_path / "shift.csv",
            ]
        )
    assert re.search(r"1 cells of .* are not in", err.value.stdout.decode("utf-8"))


def test_missing_centroid_errors(run_component, virtual_path, original_path, tmp_path):
    with pytest.raises(CalledProcessError) as err:
        run_component(
            [
                "--input",
                virtual_path,
                "--original",
                original_path,
                "--output",
                tmp_path / "shift.csv",
                "--healthy_label",
                "control",
            ]
        )
    assert re.search(r"No centroid for \['control'\]", err.value.stdout.decode("utf-8"))


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
