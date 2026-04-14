import sys
import subprocess
import re

import pytest
import mudata as mu
import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


@pytest.fixture
def source_modality():
    """Source AnnData with obs, var, obsm, varm, obsp, varp, and uns."""
    X = pd.DataFrame(
        [[1, 2, 3], [4, 5, 6]],
        index=["obs1", "obs2"],
        columns=["var1", "var2", "var3"],
    )
    obs = pd.DataFrame(
        {"cell_type": ["T-cell", "B-cell"], "score": [0.9, 0.8]},
        index=X.index,
    )
    var = pd.DataFrame(
        {"gene_name": ["GeneA", "GeneB", "GeneC"], "dispersions": [1.1, 2.2, 3.3]},
        index=X.columns,
    )
    adata = ad.AnnData(X, obs=obs, var=var)
    adata.obsm["X_pca"] = np.array([[0.1, 0.2], [0.3, 0.4]])
    adata.varm["PCs"] = np.array([[0.5, 0.6], [0.7, 0.8], [0.9, 1.0]])
    adata.obsp["distances"] = sp.csr_matrix(np.array([[0.0, 1.0], [1.0, 0.0]]))
    adata.varp["correlations"] = sp.csr_matrix(
        np.array([[1.0, 0.5, 0.2], [0.5, 1.0, 0.3], [0.2, 0.3, 1.0]])
    )
    adata.uns["method"] = "custom_pipeline"
    adata.uns["params"] = {"n_neighbors": 15}
    return adata


@pytest.fixture
def target_modality():
    """Target AnnData with the same observations and variables but different annotations."""
    X = pd.DataFrame(
        [[10, 20, 30], [40, 50, 60]],
        index=["obs1", "obs2"],
        columns=["var1", "var2", "var3"],
    )
    obs = pd.DataFrame(
        {"batch": ["A", "B"]},
        index=X.index,
    )
    var = pd.DataFrame(
        {"highly_variable": [True, False, True]},
        index=X.columns,
    )
    return ad.AnnData(X, obs=obs, var=var)


@pytest.fixture
def other_modality():
    """A second modality in the target file that should remain untouched."""
    X = pd.DataFrame(
        [[100, 200], [300, 400]],
        index=["obs1", "obs2"],
        columns=["protA", "protB"],
    )
    return ad.AnnData(X)


@pytest.fixture
def source_h5mu_path(write_mudata_to_file, source_modality):
    mdata = mu.MuData({"rna": source_modality})
    return write_mudata_to_file(mdata)


@pytest.fixture
def target_h5mu_path(write_mudata_to_file, target_modality, other_modality):
    mdata = mu.MuData({"rna": target_modality, "prot": other_modality})
    return write_mudata_to_file(mdata)


def test_move_obs_columns(
    run_component, random_path, source_h5mu_path, target_h5mu_path
):
    """Move selected .obs columns from source to target."""
    output = random_path(extension="h5mu")
    run_component(
        [
            "--input_source", str(source_h5mu_path),
            "--input_target", str(target_h5mu_path),
            "--source_modality", "rna",
            "--obs", "cell_type",
            "--output", str(output),
        ]
    )
    assert output.is_file()

    result = mu.read_h5mu(output)
    rna = result.mod["rna"]

    # Moved column is present
    assert "cell_type" in rna.obs.columns
    assert list(rna.obs["cell_type"]) == ["T-cell", "B-cell"]

    # Non-requested column was NOT moved
    assert "score" not in rna.obs.columns

    # Pre-existing target column is preserved
    assert "batch" in rna.obs.columns

    # Other modality is untouched
    assert "prot" in result.mod
    assert result.mod["prot"].n_vars == 2


def test_move_var_columns(
    run_component, random_path, source_h5mu_path, target_h5mu_path
):
    """Move selected .var columns from source to target."""
    output = random_path(extension="h5mu")
    run_component(
        [
            "--input_source", str(source_h5mu_path),
            "--input_target", str(target_h5mu_path),
            "--source_modality", "rna",
            "--var", "gene_name", "--var", "dispersions",
            "--output", str(output),
        ]
    )
    rna = mu.read_h5mu(output).mod["rna"]

    assert "gene_name" in rna.var.columns
    assert "dispersions" in rna.var.columns
    assert list(rna.var["gene_name"]) == ["GeneA", "GeneB", "GeneC"]

    # Pre-existing target column preserved
    assert "highly_variable" in rna.var.columns


def test_move_obsm_and_uns(
    run_component, random_path, source_h5mu_path, target_h5mu_path
):
    """Move .obsm and .uns keys from source to target."""
    output = random_path(extension="h5mu")
    run_component(
        [
            "--input_source", str(source_h5mu_path),
            "--input_target", str(target_h5mu_path),
            "--source_modality", "rna",
            "--obsm", "X_pca",
            "--uns", "method",
            "--output", str(output),
        ]
    )
    rna = mu.read_h5mu(output).mod["rna"]

    assert "X_pca" in rna.obsm
    np.testing.assert_array_almost_equal(
        rna.obsm["X_pca"], np.array([[0.1, 0.2], [0.3, 0.4]])
    )
    assert rna.uns["method"] == "custom_pipeline"

    # uns key not requested should not appear
    assert "params" not in rna.uns


def test_move_varm_obsp_varp(
    run_component, random_path, source_h5mu_path, target_h5mu_path
):
    """Move .varm, .obsp, and .varp keys from source to target."""
    output = random_path(extension="h5mu")
    run_component(
        [
            "--input_source", str(source_h5mu_path),
            "--input_target", str(target_h5mu_path),
            "--source_modality", "rna",
            "--varm", "PCs",
            "--obsp", "distances",
            "--varp", "correlations",
            "--output", str(output),
        ]
    )
    rna = mu.read_h5mu(output).mod["rna"]

    assert "PCs" in rna.varm
    np.testing.assert_array_almost_equal(
        rna.varm["PCs"], np.array([[0.5, 0.6], [0.7, 0.8], [0.9, 1.0]])
    )
    assert "distances" in rna.obsp
    assert rna.obsp["distances"].shape == (2, 2)
    assert "correlations" in rna.varp
    assert rna.varp["correlations"].shape == (3, 3)


def test_move_all_slots(
    run_component, random_path, source_h5mu_path, target_h5mu_path
):
    """Move at least one key from every slot type in a single invocation."""
    output = random_path(extension="h5mu")
    run_component(
        [
            "--input_source", str(source_h5mu_path),
            "--input_target", str(target_h5mu_path),
            "--source_modality", "rna",
            "--obs", "cell_type",
            "--var", "gene_name",
            "--obsm", "X_pca",
            "--varm", "PCs",
            "--obsp", "distances",
            "--varp", "correlations",
            "--uns", "method",
            "--output", str(output),
        ]
    )
    rna = mu.read_h5mu(output).mod["rna"]

    assert "cell_type" in rna.obs.columns
    assert "gene_name" in rna.var.columns
    assert "X_pca" in rna.obsm
    assert "PCs" in rna.varm
    assert "distances" in rna.obsp
    assert "correlations" in rna.varp
    assert rna.uns["method"] == "custom_pipeline"


def test_target_modality_defaults_to_source(
    run_component, random_path, source_h5mu_path, target_h5mu_path
):
    """When --target_modality is not set, it defaults to --source_modality."""
    output = random_path(extension="h5mu")
    run_component(
        [
            "--input_source", str(source_h5mu_path),
            "--input_target", str(target_h5mu_path),
            "--source_modality", "rna",
            "--obs", "cell_type",
            "--output", str(output),
        ]
    )
    result = mu.read_h5mu(output)
    assert "cell_type" in result.mod["rna"].obs.columns


def test_cross_modality_move(
    run_component, random_path, write_mudata_to_file, source_modality, other_modality
):
    """Move from rna modality in source to prot modality in target using --target_modality."""
    # Source has rna with obs columns
    source_path = write_mudata_to_file(mu.MuData({"rna": source_modality}))

    # Target has prot modality with same observations
    prot = ad.AnnData(
        pd.DataFrame(
            [[100, 200], [300, 400]],
            index=["obs1", "obs2"],
            columns=["protA", "protB"],
        )
    )
    target_path = write_mudata_to_file(mu.MuData({"prot": prot}))

    output = random_path(extension="h5mu")
    run_component(
        [
            "--input_source", str(source_path),
            "--input_target", str(target_path),
            "--source_modality", "rna",
            "--target_modality", "prot",
            "--obs", "cell_type",
            "--output", str(output),
        ]
    )
    result = mu.read_h5mu(output)
    assert "cell_type" in result.mod["prot"].obs.columns
    assert list(result.mod["prot"].obs["cell_type"]) == ["T-cell", "B-cell"]


def test_missing_obs_column_raises_error(
    run_component, random_path, source_h5mu_path, target_h5mu_path
):
    """Requesting a non-existent .obs column should fail with a clear error."""
    output = random_path(extension="h5mu")
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input_source", str(source_h5mu_path),
                "--input_target", str(target_h5mu_path),
                "--source_modality", "rna",
                "--obs", "nonexistent_column",
                "--output", str(output),
            ]
        )
    assert re.search(
        r"ValueError.*\.obs columns were not found.*nonexistent_column",
        err.value.stdout.decode("utf-8"),
    )


def test_missing_obsm_key_raises_error(
    run_component, random_path, source_h5mu_path, target_h5mu_path
):
    """Requesting a non-existent .obsm key should fail with a clear error."""
    output = random_path(extension="h5mu")
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input_source", str(source_h5mu_path),
                "--input_target", str(target_h5mu_path),
                "--source_modality", "rna",
                "--obsm", "X_umap",
                "--output", str(output),
            ]
        )
    assert re.search(
        r"ValueError.*\.obsm keys were not found.*X_umap",
        err.value.stdout.decode("utf-8"),
    )


def test_missing_source_modality_raises_error(
    run_component, random_path, source_h5mu_path, target_h5mu_path
):
    """Requesting a non-existent source modality should fail with a clear error."""
    output = random_path(extension="h5mu")
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input_source", str(source_h5mu_path),
                "--input_target", str(target_h5mu_path),
                "--source_modality", "nonexistent",
                "--obs", "cell_type",
                "--output", str(output),
            ]
        )
    assert re.search(
        r"ValueError.*nonexistent.*does not exist in source file",
        err.value.stdout.decode("utf-8"),
    )


def test_overwrites_existing_slot(
    run_component, random_path, write_mudata_to_file, source_modality
):
    """Moving a column that already exists in the target should overwrite it."""
    # Target already has a "cell_type" column with different values
    target_adata = ad.AnnData(
        pd.DataFrame(
            [[10, 20, 30], [40, 50, 60]],
            index=["obs1", "obs2"],
            columns=["var1", "var2", "var3"],
        ),
        obs=pd.DataFrame(
            {"cell_type": ["old_A", "old_B"]},
            index=["obs1", "obs2"],
        ),
    )
    source_path = write_mudata_to_file(mu.MuData({"rna": source_modality}))
    target_path = write_mudata_to_file(mu.MuData({"rna": target_adata}))

    output = random_path(extension="h5mu")
    run_component(
        [
            "--input_source", str(source_path),
            "--input_target", str(target_path),
            "--source_modality", "rna",
            "--obs", "cell_type",
            "--output", str(output),
        ]
    )
    rna = mu.read_h5mu(output).mod["rna"]
    assert list(rna.obs["cell_type"]) == ["T-cell", "B-cell"]


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
