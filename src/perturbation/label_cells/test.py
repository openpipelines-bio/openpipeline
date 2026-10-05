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
    "executable": "./target/executable/perturbation/label_cells/label_cells",
    "resources_dir": "src/utils",
    "config": "./src/perturbation/label_cells/config.vsh.yaml",
}
## VIASH END


@pytest.fixture
def input_h5mu():
    # cluster x reference grid: every combination of cluster 0-3 and
    # reference 0-2 appears exactly twice
    clusters = np.repeat(np.arange(4), 6)
    reference = np.tile(np.repeat(np.arange(3), 2), 4)
    n_obs = clusters.shape[0]
    rna = ad.AnnData(
        X=np.ones((n_obs, 2)),
        obs=pd.DataFrame(
            {
                "leiden": pd.Categorical(clusters.astype(str)),
                "reference": reference,
            },
            index=[f"cell_{i}" for i in range(n_obs)],
        ),
    )
    prot = ad.AnnData(X=np.ones((n_obs, 1)), obs=pd.DataFrame(index=rna.obs_names))
    return mu.MuData({"rna": rna, "prot": prot})


@pytest.fixture
def input_path(input_h5mu, random_h5mu_path):
    path = random_h5mu_path()
    input_h5mu.write_h5mu(path)
    return path


def test_single_clustering(run_component, input_path, random_h5mu_path):
    output_path = random_h5mu_path()
    run_component(
        [
            "--input",
            input_path,
            "--output",
            output_path,
            "--obs_cluster",
            "leiden",
            "--disease_clusters",
            "1;2",
            "--healthy_clusters",
            "3",
        ]
    )
    output = mu.read_h5mu(output_path)
    obs = output.mod["rna"].obs
    labels = obs["perturbation_group"]
    assert list(labels.cat.categories) == ["disease", "healthy", "unused"]
    assert (labels[obs["leiden"].isin(["1", "2"])] == "disease").all()
    assert (labels[obs["leiden"] == "3"] == "healthy").all()
    assert (labels[obs["leiden"] == "0"] == "unused").all()

    rules = output.mod["rna"].uns["perturbation_label_cells"]
    assert rules["n_cells"]["disease"] == 12
    assert rules["n_cells"]["healthy"] == 6
    assert rules["n_cells"]["unused"] == 6
    assert output.mod["prot"].shape == (24, 1)


def test_reference_clustering_narrows_the_rule(
    run_component, input_path, random_h5mu_path
):
    output_path = random_h5mu_path()
    run_component(
        [
            "--input",
            input_path,
            "--output",
            output_path,
            "--obs_cluster",
            "leiden",
            "--obs_reference",
            "reference",
            "--disease_clusters",
            "1",
            "--disease_reference_clusters",
            "0;2",
            "--healthy_clusters",
            "1",
            "--healthy_reference_clusters",
            "1",
            "--obs_output",
            "my_labels",
            "--output_compression",
            "gzip",
        ]
    )
    obs = mu.read_h5mu(output_path).mod["rna"].obs
    in_cluster = obs["leiden"] == "1"
    disease = in_cluster & obs["reference"].isin([0, 2])
    healthy = in_cluster & (obs["reference"] == 1)
    assert (obs.loc[disease, "my_labels"] == "disease").all()
    assert (obs.loc[healthy, "my_labels"] == "healthy").all()
    assert (obs.loc[~in_cluster, "my_labels"] == "unused").all()
    assert disease.sum() == 4 and healthy.sum() == 2


def test_overlapping_rules_error(run_component, input_path, random_h5mu_path):
    with pytest.raises(CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_path,
                "--output",
                random_h5mu_path(),
                "--obs_cluster",
                "leiden",
                "--disease_clusters",
                "1;2",
                "--healthy_clusters",
                "2",
            ]
        )
    assert re.search(r"both match 6 cells", err.value.stdout.decode("utf-8"))


def test_reference_clusters_need_reference_column(
    run_component, input_path, random_h5mu_path
):
    with pytest.raises(CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_path,
                "--output",
                random_h5mu_path(),
                "--obs_cluster",
                "leiden",
                "--disease_clusters",
                "1",
                "--disease_reference_clusters",
                "0",
                "--healthy_clusters",
                "2",
            ]
        )
    assert re.search(r"--obs_reference was not", err.value.stdout.decode("utf-8"))


def test_too_few_cells_error(run_component, input_path, random_h5mu_path):
    with pytest.raises(CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_path,
                "--output",
                random_h5mu_path(),
                "--obs_cluster",
                "leiden",
                "--disease_clusters",
                "1",
                "--healthy_clusters",
                "does_not_exist",
            ]
        )
    assert re.search(r"healthy rule matched 0 cells", err.value.stdout.decode("utf-8"))


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
