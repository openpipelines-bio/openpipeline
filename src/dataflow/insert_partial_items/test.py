import sys

import anndata as ad
import mudata
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

## VIASH START
meta = {
    "executable": "./target/executable/dataflow/insert_partial_items/insert_partial_items",
    "resources_dir": "src/utils",
    "config": "./src/dataflow/insert_partial_items/config.vsh.yaml",
}
## VIASH END


def _make_base(rng):
    obs_names = [f"cell_{i}" for i in range(20)]
    var_names = [f"gene_{i}" for i in range(15)]
    X = rng.integers(0, 100, size=(20, 15)).astype(float)
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=obs_names),
        var=pd.DataFrame(index=var_names),
    )
    adata.layers["counts"] = X.copy()
    return mudata.MuData({"rna": adata})


def _make_insert(rng):
    # Subset of base's cells and vars (base must be a superset per contract)
    obs_names = [f"cell_{i}" for i in range(15)]
    var_names = [f"gene_{i}" for i in range(12)]
    n_obs, n_vars = len(obs_names), len(var_names)

    X = rng.integers(0, 100, size=(n_obs, n_vars)).astype(float)
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(
            {
                "new_string": [f"s_{i}" for i in range(n_obs)],
                "new_integer": rng.integers(0, 100, size=n_obs),
                "new_float": rng.random(n_obs),
                "new_boolean": rng.integers(0, 2, size=n_obs).astype(bool),
            },
            index=obs_names,
        ),
        var=pd.DataFrame(
            {
                "new_string": [f"v_{i}" for i in range(n_vars)],
                "new_integer": rng.integers(0, 100, size=n_vars),
                "new_float": rng.random(n_vars),
                "new_boolean": rng.integers(0, 2, size=n_vars).astype(bool),
            },
            index=var_names,
        ),
    )
    adata.layers["logcounts"] = np.log1p(X)
    adata.obsm["X_pca"] = rng.random((n_obs, 5))
    adata.obsm["new_df"] = pd.DataFrame(
        {
            "df_col_a": rng.random(n_obs),
            "df_col_b": [f"x_{i}" for i in range(n_obs)],
        },
        index=obs_names,
    )
    adata.varm["PCs"] = rng.random((n_vars, 5))
    adata.obsp["connectivities"] = sp.random(n_obs, n_obs, density=0.1, format="csr")
    adata.obsp["distances"] = sp.random(n_obs, n_obs, density=0.1, format="csr")
    adata.uns["log1p"] = {"base": None}
    adata.uns["neighbors"] = {"params": {"n_neighbors": 5}}
    adata.uns["pca"] = {"variance_ratio": rng.random(5)}
    return mudata.MuData({"rna": adata})


@pytest.fixture
def base_mdata():
    rng = np.random.default_rng(seed=1)
    return _make_base(rng)


@pytest.fixture
def insert_mdata():
    rng = np.random.default_rng(seed=2)
    return _make_insert(rng)


@pytest.fixture
def base_path(base_mdata, random_h5mu_path):
    path = random_h5mu_path()
    base_mdata.write_h5mu(path)
    return path


@pytest.fixture
def insert_path(insert_mdata, random_h5mu_path):
    path = random_h5mu_path()
    insert_mdata.write_h5mu(path)
    return path


def get_common_indices(adata_base, adata_insert, adata_out):
    common_obs_names = adata_base.obs_names.intersection(adata_insert.obs_names)
    base_obs_idx = adata_base.obs_names.get_indexer(common_obs_names)
    insert_obs_idx = adata_insert.obs_names.get_indexer(common_obs_names)
    missing_obs_mask = ~adata_out.obs_names.isin(adata_insert.obs_names)

    common_var_names = adata_base.var_names.intersection(adata_insert.var_names)
    base_var_idx = adata_base.var_names.get_indexer(common_var_names)
    insert_var_idx = adata_insert.var_names.get_indexer(common_var_names)
    missing_var_mask = ~adata_out.var_names.isin(adata_insert.var_names)

    return (
        base_obs_idx,
        insert_obs_idx,
        missing_obs_mask,
        base_var_idx,
        insert_var_idx,
        missing_var_mask,
    )


def assert_dtype_match(output_series, insert_series):
    insert_dtype = insert_series.dtype
    output_dtype = output_series.dtype

    if pd.api.types.is_bool_dtype(insert_dtype):
        assert pd.api.types.is_bool_dtype(output_dtype)
        return

    if pd.api.types.is_integer_dtype(insert_dtype):
        assert pd.api.types.is_integer_dtype(output_dtype)
        return

    if pd.api.types.is_float_dtype(insert_dtype):
        assert pd.api.types.is_float_dtype(output_dtype)
        return

    if pd.api.types.is_string_dtype(insert_dtype) or insert_dtype is object:
        assert (
            pd.api.types.is_string_dtype(output_dtype)
            or output_dtype is object
            or isinstance(output_dtype, pd.CategoricalDtype)
        )
        return

    if isinstance(insert_dtype, pd.CategoricalDtype):
        assert isinstance(output_dtype, pd.CategoricalDtype)
        assert list(output_series.cat.categories) == list(insert_series.cat.categories)
        assert output_series.cat.ordered == insert_series.cat.ordered
        return

    assert output_dtype == insert_dtype


def assert_inserted_columns(
    output_df, insert_df, columns, base_idx, insert_idx, missing_mask
):
    for col in columns:
        assert col in output_df.columns
        out_vals = output_df[col].to_numpy()
        insert_vals = insert_df.iloc[insert_idx][col].to_numpy()
        np.testing.assert_array_equal(out_vals[base_idx], insert_vals)
        assert pd.isna(out_vals[missing_mask]).all()
        assert_dtype_match(output_df[col], insert_df[col])


def test_simple_execution(
    run_component,
    base_mdata,
    insert_mdata,
    base_path,
    insert_path,
    random_h5mu_path,
):
    output = random_h5mu_path()

    run_component(
        [
            "--input_base",
            base_path,
            "--modality_base",
            "rna",
            "--input_insert",
            insert_path,
            "--modality_insert",
            "rna",
            "--layers",
            "logcounts",
            "--obs",
            "new_string,new_integer,new_float,new_boolean",
            "--var",
            "new_string,new_integer,new_float,new_boolean",
            "--obsm",
            "X_pca,new_df",
            "--varm",
            "PCs",
            "--obsp",
            "connectivities,distances",
            "--uns",
            "log1p,neighbors,pca",
            "--output",
            output,
        ]
    )
    assert output.is_file(), "output file was not created"

    mdata_out = mudata.read_h5mu(output)
    adata_base = base_mdata["rna"]
    adata_insert = insert_mdata["rna"]
    adata_out = mdata_out["rna"]

    (
        base_obs_idx,
        insert_obs_idx,
        missing_obs_mask,
        base_var_idx,
        insert_var_idx,
        missing_var_mask,
    ) = get_common_indices(adata_base, adata_insert, adata_out)

    columns = ["new_string", "new_integer", "new_float", "new_boolean"]

    assert_inserted_columns(
        adata_out.obs,
        adata_insert.obs,
        columns,
        base_obs_idx,
        insert_obs_idx,
        missing_obs_mask,
    )

    assert_inserted_columns(
        adata_out.var,
        adata_insert.var,
        columns,
        base_var_idx,
        insert_var_idx,
        missing_var_mask,
    )

    assert "logcounts" in adata_out.layers
    out_layer = adata_out.layers["logcounts"]
    if sp.issparse(out_layer):
        out_layer = out_layer.toarray()
    insert_layer = adata_insert.layers["logcounts"][
        np.ix_(insert_obs_idx, insert_var_idx)
    ]
    if sp.issparse(insert_layer):
        insert_layer = insert_layer.toarray()
    np.testing.assert_allclose(
        out_layer[np.ix_(base_obs_idx, base_var_idx)],
        insert_layer,
        equal_nan=True,
    )
    assert np.isnan(out_layer[missing_obs_mask, :]).all()

    assert "X_pca" in adata_out.obsm
    np.testing.assert_allclose(
        adata_out.obsm["X_pca"][base_obs_idx, :],
        adata_insert.obsm["X_pca"][insert_obs_idx, :],
        equal_nan=True,
    )
    assert np.isnan(adata_out.obsm["X_pca"][missing_obs_mask, :]).all()

    assert "new_df" in adata_out.obsm
    assert isinstance(adata_out.obsm["new_df"], pd.DataFrame)
    assert list(adata_out.obsm["new_df"].index) == list(adata_out.obs_names)
    assert_inserted_columns(
        adata_out.obsm["new_df"],
        adata_insert.obsm["new_df"],
        list(adata_insert.obsm["new_df"].columns),
        base_obs_idx,
        insert_obs_idx,
        missing_obs_mask,
    )

    assert "PCs" in adata_out.varm
    np.testing.assert_allclose(
        adata_out.varm["PCs"][base_var_idx, :],
        adata_insert.varm["PCs"][insert_var_idx, :],
        equal_nan=True,
    )
    assert np.isnan(adata_out.varm["PCs"][missing_var_mask, :]).all()

    for key in ["connectivities", "distances"]:
        assert key in adata_out.obsp
        out_obsp = adata_out.obsp[key]
        insert_obsp = adata_insert.obsp[key]
        if sp.issparse(out_obsp):
            out_obsp = out_obsp.toarray()
        if sp.issparse(insert_obsp):
            insert_obsp = insert_obsp.toarray()
        np.testing.assert_allclose(
            out_obsp[np.ix_(base_obs_idx, base_obs_idx)],
            insert_obsp[np.ix_(insert_obs_idx, insert_obs_idx)],
            equal_nan=True,
        )

    for key in ["log1p", "neighbors", "pca"]:
        assert key in adata_out.uns


def test_misaligned_obs_var(
    run_component,
    base_mdata,
    insert_mdata,
    base_path,
    random_h5mu_path,
    tmp_path,
):
    """Test with misaligned obs and var between base and insert"""
    output = random_h5mu_path()

    adata_insert = insert_mdata["rna"]
    rng = np.random.default_rng(seed=3)
    adata_insert = adata_insert[
        rng.permutation(adata_insert.n_obs),
        rng.permutation(adata_insert.n_vars),
    ].copy()
    mdata_insert_reordered = mudata.MuData({"rna": adata_insert})

    insert_reordered = tmp_path / "insert_reordered.h5mu"
    mdata_insert_reordered.write_h5mu(insert_reordered)

    run_component(
        [
            "--input_base",
            base_path,
            "--modality_base",
            "rna",
            "--input_insert",
            str(insert_reordered),
            "--modality_insert",
            "rna",
            "--layers",
            "logcounts",
            "--obs",
            "new_string,new_integer,new_float,new_boolean",
            "--var",
            "new_string,new_integer,new_float,new_boolean",
            "--output",
            output,
        ]
    )
    assert output.is_file()

    mdata_out = mudata.read_h5mu(output)
    adata_out = mdata_out["rna"]
    adata_base = base_mdata["rna"]

    (
        base_obs_idx,
        insert_obs_idx,
        missing_obs_mask,
        base_var_idx,
        insert_var_idx,
        missing_var_mask,
    ) = get_common_indices(adata_base, adata_insert, adata_out)

    columns = ["new_string", "new_integer", "new_float", "new_boolean"]

    assert_inserted_columns(
        adata_out.obs,
        adata_insert.obs,
        columns,
        base_obs_idx,
        insert_obs_idx,
        missing_obs_mask,
    )

    assert_inserted_columns(
        adata_out.var,
        adata_insert.var,
        columns,
        base_var_idx,
        insert_var_idx,
        missing_var_mask,
    )

    assert "logcounts" in adata_out.layers
    out_layer = adata_out.layers["logcounts"]
    if sp.issparse(out_layer):
        out_layer = out_layer.toarray()
    insert_layer = adata_insert.layers["logcounts"][
        np.ix_(insert_obs_idx, insert_var_idx)
    ]
    if sp.issparse(insert_layer):
        insert_layer = insert_layer.toarray()
    np.testing.assert_allclose(
        out_layer[np.ix_(base_obs_idx, base_var_idx)],
        insert_layer,
        equal_nan=True,
    )


def test_auto_mode(
    run_component,
    base_mdata,
    insert_mdata,
    base_path,
    insert_path,
    random_h5mu_path,
):
    """Test the __auto__ mode to automatically include items not in base"""
    output = random_h5mu_path()

    adata_base = base_mdata["rna"]
    adata_insert = insert_mdata["rna"]

    auto_obs_cols = set(adata_insert.obs.columns) - set(adata_base.obs.columns)
    auto_var_cols = set(adata_insert.var.columns) - set(adata_base.var.columns)
    auto_layers = set(adata_insert.layers.keys()) - set(adata_base.layers.keys())
    auto_obsm = set(adata_insert.obsm.keys()) - set(adata_base.obsm.keys())
    auto_uns = set(adata_insert.uns.keys()) - set(adata_base.uns.keys())

    run_component(
        [
            "--input_base",
            base_path,
            "--modality_base",
            "rna",
            "--input_insert",
            insert_path,
            "--modality_insert",
            "rna",
            "--obs",
            "__auto__",
            "--var",
            "__auto__",
            "--layers",
            "__auto__",
            "--obsm",
            "__auto__",
            "--uns",
            "__auto__",
            "--output",
            output,
        ]
    )
    assert output.is_file(), "output file was not created"

    mdata_out = mudata.read_h5mu(output)
    adata_out = mdata_out["rna"]

    for col in auto_obs_cols:
        assert col in adata_out.obs.columns, (
            f"Auto-selected obs column '{col}' not found"
        )

    for col in auto_var_cols:
        assert col in adata_out.var.columns, (
            f"Auto-selected var column '{col}' not found"
        )

    for layer in auto_layers:
        assert layer in adata_out.layers, f"Auto-selected layer '{layer}' not found"

    for obsm_key in auto_obsm:
        assert obsm_key in adata_out.obsm, f"Auto-selected obsm '{obsm_key}' not found"

    for uns_key in auto_uns:
        assert uns_key in adata_out.uns, f"Auto-selected uns '{uns_key}' not found"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
