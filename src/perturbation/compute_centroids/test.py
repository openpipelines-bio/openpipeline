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
    "executable": "./target/executable/perturbation/compute_centroids/compute_centroids",
    "resources_dir": "src/utils",
    "config": "./src/perturbation/compute_centroids/config.vsh.yaml",
}
## VIASH END

N_OBS, N_DIMS = 30, 8


@pytest.fixture
def input_h5mu():
    rng = np.random.default_rng(2)
    groups = np.array(["disease"] * 10 + ["healthy"] * 12 + ["unused"] * 8)
    selected = np.ones(N_OBS, dtype=bool)
    selected[[0, 1, 10]] = False
    rna = ad.AnnData(
        X=np.zeros((N_OBS, 2)),
        obs=pd.DataFrame(
            {"perturbation_group": pd.Categorical(groups), "selected": selected},
            index=[f"cell_{i}" for i in range(N_OBS)],
        ),
    )
    rna.obsm["X_geneformer"] = rng.normal(size=(N_OBS, N_DIMS)).astype(np.float32)
    return mu.MuData({"rna": rna})


@pytest.fixture
def input_path(input_h5mu, random_h5mu_path):
    path = random_h5mu_path()
    input_h5mu.write_h5mu(path)
    return path


def test_median_centroids_of_filtered_groups(
    run_component, input_h5mu, input_path, random_h5mu_path
):
    output_path = random_h5mu_path()
    run_component(
        [
            "--input",
            input_path,
            "--output",
            output_path,
            "--obs_filter",
            "selected",
            "--groups",
            "disease;healthy",
        ]
    )
    centroids = mu.read_h5mu(output_path).mod["rna"].uns["perturbation_centroids"]
    assert list(centroids.index) == ["disease", "healthy"]
    assert centroids.shape == (2, N_DIMS)

    rna = input_h5mu.mod["rna"]
    for group in ["disease", "healthy"]:
        cells = (rna.obs["perturbation_group"] == group) & rna.obs["selected"]
        expected = np.median(rna.obsm["X_geneformer"][cells.to_numpy()], axis=0)
        assert np.allclose(centroids.loc[group].to_numpy(), expected, atol=1e-6)


def test_mean_of_all_groups(run_component, input_h5mu, input_path, random_h5mu_path):
    output_path = random_h5mu_path()
    run_component(
        [
            "--input",
            input_path,
            "--output",
            output_path,
            "--method",
            "mean",
            "--uns_output",
            "my_centroids",
            "--output_compression",
            "gzip",
        ]
    )
    centroids = mu.read_h5mu(output_path).mod["rna"].uns["my_centroids"]
    assert sorted(centroids.index) == ["disease", "healthy", "unused"]
    rna = input_h5mu.mod["rna"]
    cells = (rna.obs["perturbation_group"] == "unused").to_numpy()
    expected = rna.obsm["X_geneformer"][cells].astype(np.float64).mean(axis=0)
    assert np.allclose(centroids.loc["unused"].to_numpy(), expected)


def test_empty_group_errors(run_component, input_path, random_h5mu_path):
    with pytest.raises(CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_path,
                "--output",
                random_h5mu_path(),
                "--groups",
                "disease;does_not_exist",
            ]
        )
    assert re.search(r"'does_not_exist' has 0 cells", err.value.stdout.decode("utf-8"))


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
