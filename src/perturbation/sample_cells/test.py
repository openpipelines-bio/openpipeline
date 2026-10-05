import re
import sys
from subprocess import CalledProcessError

import anndata as ad
import mudata as mu
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix

## VIASH START
meta = {
    "executable": "./target/executable/perturbation/sample_cells/sample_cells",
    "resources_dir": "src/utils",
    "config": "./src/perturbation/sample_cells/config.vsh.yaml",
}
## VIASH END

GENES = [f"GENE{i}" for i in range(30)] + [
    "MT-CO1",
    "MT-ND4",
    "MT-ATP6",
    "RPS3",
    "RPL7",
    "HBB",
]
N_DISEASE, N_HEALTHY, N_UNUSED = 30, 20, 10


@pytest.fixture
def input_h5mu():
    rng = np.random.default_rng(1)
    n_obs = N_DISEASE + N_HEALTHY + N_UNUSED
    counts = rng.integers(5, 30, size=(n_obs, len(GENES)))
    # disease cells i express 10 + i // 2 of the 30 regular genes, so the
    # depth filter has something to cut
    for i in range(N_DISEASE):
        counts[i, 10 + i // 2 : 30] = 0
    # mitochondrial genes stay a small fraction
    counts[:, 30:33] = rng.integers(0, 2, size=(n_obs, 3))
    # one healthy cell is mostly mitochondrial
    counts[N_DISEASE, 30:33] = 2000

    labels = ["disease"] * N_DISEASE + ["healthy"] * N_HEALTHY + ["unused"] * N_UNUSED
    rna = ad.AnnData(
        X=csr_matrix(counts.astype(np.float32)),
        obs=pd.DataFrame(
            {"perturbation_group": pd.Categorical(labels)},
            index=[f"cell_{i:02d}" for i in range(n_obs)],
        ),
        var=pd.DataFrame(
            {"gene_symbol": GENES}, index=[f"ENSG{i:011d}" for i in range(len(GENES))]
        ),
    )
    prot = ad.AnnData(X=np.ones((n_obs, 1)), obs=pd.DataFrame(index=rna.obs_names))
    return mu.MuData({"rna": rna, "prot": prot})


@pytest.fixture
def input_path(input_h5mu, random_h5mu_path):
    path = random_h5mu_path()
    input_h5mu.write_h5mu(path)
    return path


def n_genes(adata):
    return np.asarray((adata.X > 0).sum(axis=1)).ravel()


def test_depth_filter_and_batches(run_component, input_path, random_h5mu_path):
    output_path = random_h5mu_path()
    run_component(
        [
            "--input",
            input_path,
            "--output",
            output_path,
            "--var_gene_names",
            "gene_symbol",
            "--skip_qc",
            "--min_genes",
            "25",
            "--batch_size",
            "3",
        ]
    )
    rna = mu.read_h5mu(output_path).mod["rna"]
    obs = rna.obs

    is_disease = (obs["perturbation_group"] == "disease").to_numpy()
    is_healthy = (obs["perturbation_group"] == "healthy").to_numpy()
    deep = n_genes(rna) > 25
    expected_perturb = is_disease & deep
    assert expected_perturb.sum() > 3, "fixture should keep several disease cells"
    assert np.array_equal(obs["perturbation_perturb"].to_numpy(), expected_perturb)
    assert np.array_equal(
        obs["perturbation_selected"].to_numpy(), expected_perturb | is_healthy
    )

    batches = obs.loc[expected_perturb, "perturbation_batch"].astype(str)
    sizes = batches.value_counts().sort_index()
    assert sizes.iloc[:-1].eq(3).all() and sizes.iloc[-1] <= 3
    assert sizes.sum() == expected_perturb.sum()
    assert list(batches.iloc[:3]) == ["batch_0001"] * 3
    assert (obs.loc[~expected_perturb, "perturbation_batch"].astype(str) == "").all()


def test_sampling_is_seeded(run_component, input_path, random_h5mu_path):
    def sample(seed):
        output_path = random_h5mu_path()
        run_component(
            [
                "--input",
                input_path,
                "--output",
                output_path,
                "--skip_qc",
                "--min_genes",
                "0",
                "--n_disease_cells",
                "8",
                "--seed",
                str(seed),
            ]
        )
        obs = mu.read_h5mu(output_path).mod["rna"].obs
        return obs.index[obs["perturbation_perturb"]].tolist()

    first = sample(45)
    assert len(first) == 8
    assert first == sample(45)
    assert first != sample(7)


def test_qc_drops_mitochondrial_cell(run_component, input_path, random_h5mu_path):
    output_path = random_h5mu_path()
    run_component(
        [
            "--input",
            input_path,
            "--output",
            output_path,
            "--var_gene_names",
            "gene_symbol",
            "--min_genes",
            "0",
            "--output_compression",
            "gzip",
        ]
    )
    obs = mu.read_h5mu(output_path).mod["rna"].obs
    mt_cell = f"cell_{N_DISEASE:02d}"
    assert not obs.loc[mt_cell, "perturbation_selected"]
    healthy = obs["perturbation_group"] == "healthy"
    assert obs.loc[healthy, "perturbation_selected"].sum() >= N_HEALTHY - 5


def test_no_mitochondrial_genes_error(run_component, input_path, random_h5mu_path):
    # the .var index holds Ensembl ids, so no gene starts with MT-
    with pytest.raises(CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_path,
                "--output",
                random_h5mu_path(),
                "--min_genes",
                "0",
            ]
        )
    assert re.search(
        r"No mitochondrial gene \('MT-' prefix\) found in the \.var index",
        err.value.stdout.decode("utf-8"),
    )


def test_too_few_disease_cells_error(run_component, input_path, random_h5mu_path):
    with pytest.raises(CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_path,
                "--output",
                random_h5mu_path(),
                "--skip_qc",
                "--min_genes",
                "1000",
            ]
        )
    assert re.search(
        r"Only 0 disease cells are selected", err.value.stdout.decode("utf-8")
    )


def test_missing_labels_error(run_component, input_path, random_h5mu_path):
    with pytest.raises(CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_path,
                "--output",
                random_h5mu_path(),
                "--obs_group",
                "not_labelled",
            ]
        )
    assert re.search(
        r"perturbation/label_cells first", err.value.stdout.decode("utf-8")
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
