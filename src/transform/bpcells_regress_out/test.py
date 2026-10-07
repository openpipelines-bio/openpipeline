import sys
import subprocess
import pytest
import h5py
import mudata as mu
import numpy as np
from scipy.sparse import issparse

## VIASH START
meta = {"name": "bpcells_regress_out", "resources_dir": "resources_test/"}
## VIASH END


@pytest.fixture
def input_h5mu_path():
    return f"{meta['resources_dir']}/pbmc_1k_protein_v3_mms.h5mu"


@pytest.fixture
def output_h5mu_path(tmp_path):
    return tmp_path / "output.h5mu"


def test_regress_out(run_component, input_h5mu_path, output_h5mu_path):
    # execute command
    cmd_pars = [
        "--input",
        input_h5mu_path,
        "--output",
        output_h5mu_path,
        "--obs_keys",
        "total_counts",
    ]
    run_component(cmd_pars)

    assert output_h5mu_path.is_file(), "No output was created."

    mu_input = mu.read_h5mu(input_h5mu_path)
    mu_output = mu.read_h5mu(output_h5mu_path)

    assert "rna" in mu_output.mod, 'Output should contain data.mod["rna"].'
    assert "prot" in mu_output.mod, 'Output should contain data.mod["prot"].'

    rna_in = mu_input.mod["rna"]
    rna_out = mu_output.mod["rna"]
    prot_in = mu_input.mod["prot"]
    prot_out = mu_output.mod["prot"]

    assert rna_in.shape == rna_out.shape, "Should have same shape as before"
    assert prot_in.shape == prot_out.shape, "Should have same shape as before"

    assert np.mean(rna_in.X) != np.mean(rna_out.layers["regressed"]), (
        "RNA expression should have changed"
    )
    assert np.mean(prot_in.X) == np.mean(prot_out.X), (
        "Protein expression should remain the same"
    )


def test_regress_out_output_compression(
    run_component, input_h5mu_path, output_h5mu_path
):
    cmd_pars = [
        "--input",
        input_h5mu_path,
        "--output",
        output_h5mu_path,
        "--obs_keys",
        "total_counts",
        "--output_layer_compression",
        "4",
    ]
    run_component(cmd_pars)

    with h5py.File(output_h5mu_path, "r") as h5:
        data = h5["mod/rna/layers/regressed/data"]
        assert data.compression == "gzip", "Output layer should be gzip compressed"
        assert data.compression_opts == 4, "Output layer should use gzip level 4"


def test_regress_out_with_layers(run_component, input_h5mu_path, output_h5mu_path):
    # execute command
    cmd_pars = [
        "--input",
        input_h5mu_path,
        "--output",
        output_h5mu_path,
        "--obs_keys",
        "total_counts",
        "--input_layer",
        "log_normalized",
    ]
    run_component(cmd_pars)

    rna_in = mu.read_h5ad(input_h5mu_path, mod="rna")
    rna_out = mu.read_h5ad(output_h5mu_path, mod="rna")

    assert np.mean(rna_in.layers["log_normalized"]) != np.mean(
        rna_out.layers["regressed"]
    ), "RNA expression should have changed"


def test_regress_out_hvg(run_component, input_h5mu_path, output_h5mu_path, tmp_path):
    base_pars = [
        "--input",
        input_h5mu_path,
        "--obs_keys",
        "total_counts",
        "--input_layer",
        "log_normalized",
    ]
    run_component(
        base_pars + ["--output", output_h5mu_path, "--var_input", "filter_with_hvg"]
    )
    all_genes_path = tmp_path / "all_genes.h5mu"
    run_component(base_pars + ["--output", all_genes_path])

    rna_in = mu.read_h5mu(input_h5mu_path).mod["rna"]
    rna_out = mu.read_h5mu(output_h5mu_path).mod["rna"]
    rna_all_genes = mu.read_h5mu(all_genes_path).mod["rna"]
    hvg = rna_in.var["filter_with_hvg"].to_numpy()

    assert rna_in.shape == rna_out.shape, "Should have same shape as before"

    output_layer = rna_out.layers["regressed"]
    assert issparse(output_layer), "Output layer should be a sparse matrix"
    output_matrix = output_layer.toarray()

    assert not np.any(output_matrix[:, ~hvg]), "Non-selected genes should be set to 0"
    assert output_layer.nnz <= rna_in.n_obs * hvg.sum(), (
        "Only values of selected genes should be stored"
    )
    np.testing.assert_allclose(
        output_matrix[:, hvg],
        rna_all_genes.layers["regressed"].toarray()[:, hvg],
        err_msg="Selected genes should be regressed as when using all genes",
    )


def test_regress_out_existing_output_layer(
    run_component, input_h5mu_path, output_h5mu_path
):
    cmd_pars = [
        "--input",
        input_h5mu_path,
        "--output",
        output_h5mu_path,
        "--obs_keys",
        "total_counts",
        "--output_layer",
        "log_normalized",
    ]
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(cmd_pars)
    assert "Output layer log_normalized already exists" in err.value.stdout.decode(
        "utf-8"
    )


def reference_pca(regressed, num_components, max_value=None):
    # Same steps as sc.pp.scale(max_value=...) followed by sc.tl.pca
    std = regressed.std(axis=0, ddof=1)
    std[std == 0] = 1
    scaled = (regressed - regressed.mean(axis=0)) / std
    if max_value is not None:
        scaled = np.clip(scaled, -max_value, max_value)
    centered = scaled - scaled.mean(axis=0)
    u, s, vt = np.linalg.svd(centered, full_matrices=False)
    variance = s**2 / (centered.shape[0] - 1)
    total_variance = centered.var(axis=0, ddof=1).sum()
    return (
        u[:, :num_components] * s[:num_components],
        vt[:num_components].T,
        variance[:num_components],
        variance[:num_components] / total_variance,
    )


def assert_equal_up_to_sign(actual, expected, err_msg):
    # Components are compared by their relative error: the iterative SVD solver
    # converges less tightly elementwise on the last component when its
    # singular value is close to the next one
    signs = np.sign(np.sum(actual * expected, axis=0))
    relative_error = np.linalg.norm(actual * signs - expected, axis=0) / np.linalg.norm(
        expected, axis=0
    )
    assert np.all(relative_error < 1e-5), f"{err_msg}: {relative_error}"


@pytest.mark.parametrize("max_value", [None, 2.0])
def test_regress_out_pca(run_component, input_h5mu_path, output_h5mu_path, max_value):
    cmd_pars = [
        "--input",
        input_h5mu_path,
        "--output",
        output_h5mu_path,
        "--obs_keys",
        "total_counts",
        "--input_layer",
        "log_normalized",
        "--var_input",
        "filter_with_hvg",
        "--obsm_pca_output",
        "regressed_pca",
        "--num_components",
        "5",
    ]
    if max_value is not None:
        cmd_pars += ["--scale_max_value", str(max_value)]
    run_component(cmd_pars)

    rna_out = mu.read_h5ad(output_h5mu_path, mod="rna")
    hvg = rna_out.var["filter_with_hvg"].to_numpy()

    embedding = rna_out.obsm["regressed_pca"]
    loadings = rna_out.varm["regressed_pca_loadings"]
    pca_variance = rna_out.uns["regressed_pca_variance"]
    assert embedding.shape == (rna_out.n_obs, 5)
    assert loadings.shape == (rna_out.n_vars, 5)
    assert not np.any(loadings[~hvg]), "Non-selected genes should have zero loadings"

    exp_embedding, exp_loadings, exp_variance, exp_variance_ratio = reference_pca(
        rna_out.layers["regressed"].toarray()[:, hvg], 5, max_value
    )
    assert_equal_up_to_sign(embedding, exp_embedding, "Embedding should match")
    assert_equal_up_to_sign(loadings[hvg], exp_loadings, "Loadings should match")
    np.testing.assert_allclose(pca_variance["variance"], exp_variance, rtol=1e-5)
    np.testing.assert_allclose(
        pca_variance["variance_ratio"], exp_variance_ratio, rtol=1e-5
    )


def test_regress_out_pca_overwrite(
    run_component, input_h5mu_path, output_h5mu_path, tmp_path
):
    base_pars = [
        "--input",
        input_h5mu_path,
        "--obs_keys",
        "total_counts",
        "--input_layer",
        "log_normalized",
        "--var_input",
        "filter_with_hvg",
    ]
    run_component(
        base_pars
        + [
            "--output",
            output_h5mu_path,
            "--obsm_pca_output",
            "X_pca",
            "--varm_pca_output",
            "pca_loadings",
            "--uns_pca_output",
            "pca_variance",
            "--num_components",
            "5",
            "--overwrite",
        ]
    )
    regressed_only_path = tmp_path / "regressed_only.h5mu"
    run_component(base_pars + ["--output", regressed_only_path])

    rna_out = mu.read_h5ad(output_h5mu_path, mod="rna")
    rna_regressed_only = mu.read_h5ad(regressed_only_path, mod="rna")

    assert rna_out.obsm["X_pca"].shape == (
        rna_out.n_obs,
        5,
    ), "Existing .obsm slot should be overwritten"
    assert rna_out.varm["pca_loadings"].shape == (
        rna_out.n_vars,
        5,
    ), "Existing .varm slot should be overwritten"
    assert len(rna_out.uns["pca_variance"]["variance"]) == 5, (
        "Existing .uns slot should be overwritten"
    )
    np.testing.assert_allclose(
        rna_out.layers["regressed"].toarray(),
        rna_regressed_only.layers["regressed"].toarray(),
        err_msg="Output layer should not be affected by running the PCA",
    )


def test_regress_out_pca_existing_slot(
    run_component, input_h5mu_path, output_h5mu_path
):
    cmd_pars = [
        "--input",
        input_h5mu_path,
        "--output",
        output_h5mu_path,
        "--obs_keys",
        "total_counts",
        "--obsm_pca_output",
        "X_pca",
    ]
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(cmd_pars)
    assert "mod/rna/obsm/X_pca already exist" in err.value.stdout.decode("utf-8")
    assert "--overwrite" in err.value.stdout.decode("utf-8")


def test_regress_out_pca_without_output_layer(
    run_component, input_h5mu_path, output_h5mu_path
):
    cmd_pars = [
        "--input",
        input_h5mu_path,
        "--output",
        output_h5mu_path,
        "--obs_keys",
        "total_counts",
        "--output_layer",
        "",
        "--obsm_pca_output",
        "regressed_pca",
        "--num_components",
        "5",
    ]
    run_component(cmd_pars)

    rna_in = mu.read_h5ad(input_h5mu_path, mod="rna")
    rna_out = mu.read_h5ad(output_h5mu_path, mod="rna")
    assert set(rna_out.layers.keys()) == set(rna_in.layers.keys()), (
        "No output layer should be written when --output_layer is empty"
    )
    assert rna_out.obsm["regressed_pca"].shape == (rna_out.n_obs, 5)


def test_regress_out_no_output_requested(
    run_component, input_h5mu_path, output_h5mu_path
):
    cmd_pars = [
        "--input",
        input_h5mu_path,
        "--output",
        output_h5mu_path,
        "--obs_keys",
        "total_counts",
        "--output_layer",
        "",
    ]
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(cmd_pars)
    assert (
        "At least one of --output_layer or --obsm_pca_output must be provided"
        in err.value.stdout.decode("utf-8")
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
