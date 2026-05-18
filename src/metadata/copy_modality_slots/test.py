import re
import sys
import uuid
from subprocess import CalledProcessError

import anndata as ad
import mudata
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

## VIASH START
meta = {
    "executable": "./target/executable/metadata/copy_modality_slots/copy_modality_slots",
    "resources_dir": "src/utils",
    "config": "./src/metadata/copy_modality_slots/config.vsh.yaml",
}
## VIASH END


def _make_adata(obs_names, var_names, rng, *, obs_extra=None, var_extra=None):
    n_obs, n_vars = len(obs_names), len(var_names)
    X = rng.integers(0, 100, size=(n_obs, n_vars)).astype(float)
    obs = pd.DataFrame(obs_extra or {}, index=pd.Index(obs_names))
    var = pd.DataFrame(var_extra or {}, index=pd.Index(var_names))
    return ad.AnnData(X=X, obs=obs, var=var)


@pytest.fixture
def rng():
    return np.random.default_rng(seed=1)


@pytest.fixture
def target_mdata(rng):
    obs_names = [f"cell_{i}" for i in range(10)]
    var_names = [f"gene_{i}" for i in range(8)]
    return mudata.MuData({"rna": _make_adata(obs_names, var_names, rng)})


@pytest.fixture
def source_mdata_equal(rng):
    # Same index sets as target, but different obs/var content.
    obs_names = [f"cell_{i}" for i in range(10)]
    var_names = [f"gene_{i}" for i in range(8)]
    obs_extra = {
        "label": pd.Categorical(
            ["a", "b", "a", "c", "b", "a", "c", "b", "a", "c"],
            categories=["a", "b", "c"],
        ),
        "score": rng.random(10),
        "count": rng.integers(0, 100, size=10),
        "is_doublet": rng.integers(0, 2, size=10).astype(bool),
        "name": [f"n_{i}" for i in range(10)],
    }
    var_extra = {
        "highly_variable": rng.integers(0, 2, size=8).astype(bool),
        "mean": rng.random(8),
    }
    adata = _make_adata(
        obs_names, var_names, rng, obs_extra=obs_extra, var_extra=var_extra
    )
    adata.layers["logcounts"] = np.log1p(adata.X)
    adata.layers["sparse_layer"] = sp.random(
        10, 8, density=0.2, format="csr", random_state=1
    )
    adata.obsm["X_pca"] = rng.random((10, 5))
    adata.obsm["sparse_obsm"] = sp.random(
        10, 4, density=0.2, format="csr", random_state=2
    )
    adata.obsm["df_obsm"] = pd.DataFrame(
        {"a": rng.random(10), "b": [f"x_{i}" for i in range(10)]},
        index=pd.Index(obs_names),
    )
    adata.varm["PCs"] = rng.random((8, 5))
    adata.obsp["connectivities"] = sp.random(
        10, 10, density=0.1, format="csr", random_state=3
    )
    adata.varp["correlations"] = rng.random((8, 8))
    adata.uns["pca"] = {"variance_ratio": rng.random(5).tolist()}
    adata.uns["params"] = {"n_neighbors": 5}
    return mudata.MuData({"rna": adata})


@pytest.fixture
def source_mdata_subset(rng):
    # Strict subset of target (cells/genes 0..6 / 0..5 of 10/8).
    obs_names = [f"cell_{i}" for i in range(7)]
    var_names = [f"gene_{i}" for i in range(6)]
    obs_extra = {
        "label": pd.Categorical(
            ["a", "b", "c", "a", "b", "c", "a"], categories=["a", "b", "c"]
        ),
        "score": rng.random(7),
    }
    var_extra = {"highly_variable": rng.integers(0, 2, size=6).astype(bool)}
    adata = _make_adata(
        obs_names, var_names, rng, obs_extra=obs_extra, var_extra=var_extra
    )
    adata.layers["logcounts"] = np.log1p(adata.X)
    adata.obsm["X_pca"] = rng.random((7, 5))
    adata.obsp["connectivities"] = sp.random(
        7, 7, density=0.2, format="csr", random_state=4
    )
    adata.uns["pca"] = {"variance_ratio": rng.random(5).tolist()}
    return mudata.MuData({"rna": adata})


@pytest.fixture
def write_mdata(tmp_path):
    def _write(mdata):
        path = tmp_path / f"{uuid.uuid4()}.h5mu"
        mdata.write_h5mu(path)
        return path

    return _write


# ----------------------------------------------------------------------------
# Strict mode: source and target have equal obs/var index sets.
# ----------------------------------------------------------------------------


def test_strict_copy_all_slot_kinds(
    run_component, target_mdata, source_mdata_equal, write_mdata, tmp_path
):
    output = tmp_path / "out.h5mu"
    source_path = write_mdata(source_mdata_equal)
    target_path = write_mdata(target_mdata)

    run_component(
        [
            "--input_source",
            source_path,
            "--input_target",
            target_path,
            "--obs",
            "label,score,count,is_doublet,name",
            "--var",
            "highly_variable,mean",
            "--layers",
            "logcounts,sparse_layer",
            "--obsm",
            "X_pca,sparse_obsm,df_obsm",
            "--varm",
            "PCs",
            "--obsp",
            "connectivities",
            "--varp",
            "correlations",
            "--uns",
            "pca,params",
            "--output",
            output,
        ]
    )

    out = mudata.read_h5mu(output)["rna"]
    src = source_mdata_equal["rna"]

    # .obs / .var dtypes preserved appropriately
    assert out.obs["label"].dtype == "category"
    assert list(out.obs["label"].cat.categories) == ["a", "b", "c"]
    np.testing.assert_array_equal(
        out.obs["label"].astype(str).to_numpy(), src.obs["label"].astype(str).to_numpy()
    )
    np.testing.assert_array_equal(
        out.obs["count"].to_numpy(), src.obs["count"].to_numpy()
    )
    np.testing.assert_array_equal(
        out.obs["is_doublet"].to_numpy(), src.obs["is_doublet"].to_numpy()
    )
    np.testing.assert_array_equal(
        out.obs["name"].astype(str).to_numpy(), src.obs["name"].astype(str).to_numpy()
    )
    np.testing.assert_allclose(out.obs["score"].to_numpy(), src.obs["score"].to_numpy())
    np.testing.assert_array_equal(
        out.var["highly_variable"].to_numpy(), src.var["highly_variable"].to_numpy()
    )

    # .layers - dense and sparse
    np.testing.assert_allclose(out.layers["logcounts"], src.layers["logcounts"])
    assert sp.issparse(out.layers["sparse_layer"])
    np.testing.assert_allclose(
        out.layers["sparse_layer"].toarray(), src.layers["sparse_layer"].toarray()
    )

    # .obsm - dense, sparse, DataFrame
    np.testing.assert_allclose(out.obsm["X_pca"], src.obsm["X_pca"])
    assert sp.issparse(out.obsm["sparse_obsm"])
    np.testing.assert_allclose(
        out.obsm["sparse_obsm"].toarray(), src.obsm["sparse_obsm"].toarray()
    )
    assert isinstance(out.obsm["df_obsm"], pd.DataFrame)

    # .varm / .obsp / .varp
    np.testing.assert_allclose(out.varm["PCs"], src.varm["PCs"])
    assert sp.issparse(out.obsp["connectivities"])
    np.testing.assert_allclose(
        out.obsp["connectivities"].toarray(), src.obsp["connectivities"].toarray()
    )
    np.testing.assert_allclose(out.varp["correlations"], src.varp["correlations"])

    # .uns
    np.testing.assert_allclose(
        np.asarray(out.uns["pca"]["variance_ratio"]),
        np.asarray(src.uns["pca"]["variance_ratio"]),
    )
    assert out.uns["params"]["n_neighbors"] == 5


def test_strict_reorder_source(
    run_component, target_mdata, source_mdata_equal, write_mdata, tmp_path, rng
):
    # Same index sets, but source rows shuffled - values must end up correctly
    # aligned to target order in the output.
    src = source_mdata_equal["rna"]
    shuffled = src[rng.permutation(src.n_obs), rng.permutation(src.n_vars)].copy()
    source_shuffled = mudata.MuData({"rna": shuffled})

    output = tmp_path / "out.h5mu"
    run_component(
        [
            "--input_source",
            write_mdata(source_shuffled),
            "--input_target",
            write_mdata(target_mdata),
            "--obs",
            "score",
            "--obsm",
            "X_pca",
            "--output",
            output,
        ]
    )

    out = mudata.read_h5mu(output)["rna"]
    # Values should match the original (non-shuffled) ones, since the target
    # order is canonical.
    for i, name in enumerate(out.obs_names):
        assert out.obs["score"].iloc[i] == src.obs.loc[name, "score"]
    for i, name in enumerate(out.obs_names):
        np.testing.assert_allclose(
            out.obsm["X_pca"][i], src.obsm["X_pca"][src.obs_names.get_loc(name)]
        )


def test_strict_errors_when_source_smaller(
    run_component, target_mdata, source_mdata_subset, write_mdata, tmp_path
):
    with pytest.raises(CalledProcessError) as err:
        run_component(
            [
                "--input_source",
                write_mdata(source_mdata_subset),
                "--input_target",
                write_mdata(target_mdata),
                "--obs",
                "label",
                "--output",
                tmp_path / "out.h5mu",
            ]
        )
    assert re.search(
        r"differ in size.*allow_partial", err.value.stdout.decode("utf-8"), re.S
    )


def test_errors_when_source_has_extra_cells(
    run_component, target_mdata, source_mdata_equal, write_mdata, tmp_path, rng
):
    # Add a cell to source that's not in target - should error in both strict
    # and partial modes (source must be a subset of target).
    src = source_mdata_equal["rna"]
    extra = src[:1].copy()
    extra.obs_names = pd.Index(["extra_cell"])
    src2 = ad.concat([src, extra], axis=0)
    source_mod = mudata.MuData({"rna": src2})
    with pytest.raises(CalledProcessError) as err:
        run_component(
            [
                "--input_source",
                write_mdata(source_mod),
                "--input_target",
                write_mdata(target_mdata),
                "--obs",
                "score",
                "--allow_partial",
                "--output",
                tmp_path / "out.h5mu",
            ]
        )
    assert re.search(r"not present in target", err.value.stdout.decode("utf-8"))


# ----------------------------------------------------------------------------
# Partial mode: source is a strict subset of target.
# ----------------------------------------------------------------------------


def test_partial_subset_nan_fills(
    run_component, target_mdata, source_mdata_subset, write_mdata, tmp_path
):
    output = tmp_path / "out.h5mu"
    src = source_mdata_subset["rna"]
    tgt = target_mdata["rna"]

    run_component(
        [
            "--input_source",
            write_mdata(source_mdata_subset),
            "--input_target",
            write_mdata(target_mdata),
            "--obs",
            "label,score",
            "--var",
            "highly_variable",
            "--layers",
            "logcounts",
            "--obsm",
            "X_pca",
            "--obsp",
            "connectivities",
            "--uns",
            "pca",
            "--allow_partial",
            "--output",
            output,
        ]
    )

    out = mudata.read_h5mu(output)["rna"]

    # obs columns: present rows match source, missing rows are NaN/<NA>.
    present_mask = tgt.obs_names.isin(src.obs_names)
    for i, name in enumerate(out.obs_names):
        if present_mask[i]:
            assert out.obs["score"].iloc[i] == src.obs.loc[name, "score"]
        else:
            assert pd.isna(out.obs["score"].iloc[i])
            assert pd.isna(out.obs["label"].iloc[i])

    # var: missing genes get NaN
    var_present_mask = tgt.var_names.isin(src.var_names)
    for i, name in enumerate(out.var_names):
        if var_present_mask[i]:
            assert (
                out.var["highly_variable"].iloc[i]
                == src.var.loc[name, "highly_variable"]
            )
        else:
            assert pd.isna(out.var["highly_variable"].iloc[i])

    # layers: present block matches source, rest NaN
    out_layer = out.layers["logcounts"]
    assert out_layer.shape == (tgt.n_obs, tgt.n_vars)
    for i, name in enumerate(out.obs_names):
        for j, vname in enumerate(out.var_names):
            if name in src.obs_names and vname in src.var_names:
                expected = src.layers["logcounts"][
                    src.obs_names.get_loc(name), src.var_names.get_loc(vname)
                ]
                np.testing.assert_allclose(out_layer[i, j], expected)
            else:
                assert np.isnan(out_layer[i, j])

    # obsm dense: missing rows are NaN
    assert out.obsm["X_pca"].shape == (tgt.n_obs, 5)
    for i, name in enumerate(out.obs_names):
        if name in src.obs_names:
            np.testing.assert_allclose(
                out.obsm["X_pca"][i],
                src.obsm["X_pca"][src.obs_names.get_loc(name)],
            )
        else:
            assert np.isnan(out.obsm["X_pca"][i]).all()

    # obsp sparse: dimensions expanded
    assert out.obsp["connectivities"].shape == (tgt.n_obs, tgt.n_obs)

    # uns is just copied
    np.testing.assert_allclose(
        np.asarray(out.uns["pca"]["variance_ratio"]),
        np.asarray(src.uns["pca"]["variance_ratio"]),
    )


# ----------------------------------------------------------------------------
# __auto__ token
# ----------------------------------------------------------------------------


def test_auto_token_only_copies_new_keys(
    run_component, target_mdata, source_mdata_equal, write_mdata, tmp_path
):
    # Plant an existing key in the target's obs and obsm to confirm __auto__
    # skips it.
    tgt = target_mdata["rna"]
    tgt.obs["label"] = "pre_existing"  # exists in source as Categorical
    tgt.obsm["X_pca"] = np.zeros((tgt.n_obs, 5))

    output = tmp_path / "out.h5mu"
    run_component(
        [
            "--input_source",
            write_mdata(source_mdata_equal),
            "--input_target",
            write_mdata(target_mdata),
            "--obs",
            "__auto__",
            "--obsm",
            "__auto__",
            "--uns",
            "__auto__",
            "--output",
            output,
        ]
    )

    out = mudata.read_h5mu(output)["rna"]
    # 'label' was already in target -> not copied; pre_existing values remain.
    assert (out.obs["label"] == "pre_existing").all()
    # 'score' was not in target -> copied.
    assert "score" in out.obs.columns
    # X_pca was already in target -> not copied; remains zeros.
    np.testing.assert_array_equal(out.obsm["X_pca"], np.zeros((tgt.n_obs, 5)))
    # 'pca' and 'params' uns keys were absent -> copied.
    assert "pca" in out.uns and "params" in out.uns


# ----------------------------------------------------------------------------
# Overwrite policy
# ----------------------------------------------------------------------------


def test_overwrite_errors_by_default(
    run_component, target_mdata, source_mdata_equal, write_mdata, tmp_path
):
    target_mdata["rna"].obs["score"] = 0.0  # collides with source
    with pytest.raises(CalledProcessError) as err:
        run_component(
            [
                "--input_source",
                write_mdata(source_mdata_equal),
                "--input_target",
                write_mdata(target_mdata),
                "--obs",
                "score",
                "--output",
                tmp_path / "out.h5mu",
            ]
        )
    assert re.search(
        r"already exist.*allow_overwrite", err.value.stdout.decode("utf-8"), re.S
    )


def test_allow_overwrite_warns_and_overwrites(
    run_component, target_mdata, source_mdata_equal, write_mdata, tmp_path
):
    target_mdata["rna"].obs["score"] = 0.0
    output = tmp_path / "out.h5mu"
    run_component(
        [
            "--input_source",
            write_mdata(source_mdata_equal),
            "--input_target",
            write_mdata(target_mdata),
            "--obs",
            "score",
            "--allow_overwrite",
            "--output",
            output,
        ]
    )
    out = mudata.read_h5mu(output)["rna"]
    src = source_mdata_equal["rna"]
    np.testing.assert_allclose(out.obs["score"].to_numpy(), src.obs["score"].to_numpy())


# ----------------------------------------------------------------------------
# --var_match_column
# ----------------------------------------------------------------------------


def test_var_match_column_rewrites_source_var_index(
    run_component, target_mdata, source_mdata_equal, write_mdata, tmp_path
):
    # Sanitise the source's var index, but preserve the original name in a
    # column the component will use to match against the target.
    src = source_mdata_equal["rna"]
    original = list(src.var_names)
    src.var["_ori_var_index"] = original
    src.var_names = pd.Index([f"sanitised_{i}" for i in range(src.n_vars)])

    output = tmp_path / "out.h5mu"
    run_component(
        [
            "--input_source",
            write_mdata(source_mdata_equal),
            "--input_target",
            write_mdata(target_mdata),
            "--var_match_column",
            "_ori_var_index",
            "--var",
            "highly_variable",
            "--output",
            output,
        ]
    )

    out = mudata.read_h5mu(output)["rna"]
    # Values from source should now be aligned to target's gene names.
    for vname in out.var_names:
        expected = src.var.loc[
            src.var.index[src.var["_ori_var_index"] == vname][0], "highly_variable"
        ]
        assert out.var.loc[vname, "highly_variable"] == expected


# ----------------------------------------------------------------------------
# Cross-modality
# ----------------------------------------------------------------------------


def test_cross_modality_copy(
    run_component, source_mdata_equal, write_mdata, tmp_path, rng
):
    # Build a target whose 'prot' modality has the same cells/genes as
    # source's 'rna' modality.
    obs_names = list(source_mdata_equal["rna"].obs_names)
    var_names = list(source_mdata_equal["rna"].var_names)
    prot = _make_adata(obs_names, var_names, rng)
    target = mudata.MuData({"prot": prot})

    output = tmp_path / "out.h5mu"
    run_component(
        [
            "--input_source",
            write_mdata(source_mdata_equal),
            "--source_modality",
            "rna",
            "--input_target",
            write_mdata(target),
            "--target_modality",
            "prot",
            "--obs",
            "score",
            "--output",
            output,
        ]
    )

    out = mudata.read_h5mu(output)
    assert "prot" in out.mod
    np.testing.assert_allclose(
        out["prot"].obs["score"].to_numpy(),
        source_mdata_equal["rna"].obs["score"].to_numpy(),
    )


# ----------------------------------------------------------------------------
# Error paths
# ----------------------------------------------------------------------------


def test_missing_modality_errors(
    run_component, source_mdata_equal, write_mdata, target_mdata, tmp_path
):
    with pytest.raises(CalledProcessError) as err:
        run_component(
            [
                "--input_source",
                write_mdata(source_mdata_equal),
                "--source_modality",
                "doesnotexist",
                "--input_target",
                write_mdata(target_mdata),
                "--obs",
                "score",
                "--output",
                tmp_path / "out.h5mu",
            ]
        )
    assert re.search(r"doesnotexist.*does not exist", err.value.stdout.decode("utf-8"))


def test_missing_key_in_source_errors(
    run_component, source_mdata_equal, write_mdata, target_mdata, tmp_path
):
    with pytest.raises(CalledProcessError) as err:
        run_component(
            [
                "--input_source",
                write_mdata(source_mdata_equal),
                "--input_target",
                write_mdata(target_mdata),
                "--obs",
                "no_such_column",
                "--output",
                tmp_path / "out.h5mu",
            ]
        )
    output = err.value.stdout.decode("utf-8")
    assert "not found in source" in output
    assert "no_such_column" in output


def test_compression(
    run_component, target_mdata, source_mdata_equal, write_mdata, tmp_path
):
    output = tmp_path / "out.h5mu"
    run_component(
        [
            "--input_source",
            write_mdata(source_mdata_equal),
            "--input_target",
            write_mdata(target_mdata),
            "--obs",
            "score",
            "--output_compression",
            "gzip",
            "--output",
            output,
        ]
    )
    assert output.is_file()
    out = mudata.read_h5mu(output)["rna"]
    np.testing.assert_allclose(
        out.obs["score"].to_numpy(),
        source_mdata_equal["rna"].obs["score"].to_numpy(),
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
