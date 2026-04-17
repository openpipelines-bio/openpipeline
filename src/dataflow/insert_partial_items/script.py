import sys
from warnings import warn

import mudata
import numpy as np
import pandas as pd
import scipy.sparse as sp

################################################################################
# VIASH
################################################################################

## VIASH START
par = {
    "input_base": "input_base.zarr",
    "modality_base": "RNA",
    "input_insert": "input_insert.zarr",
    "modality_insert": "RNA",
    "output": "output.zarr",
    "obs": ["example_col"],
    "var": ["example_var"],
    "layers": ["example_layer"],
    "obsm": ["example_obsm"],
    "varm": ["example_varm"],
    "obsp": ["example_obsp"],
    "varp": ["example_varp"],
    "uns": ["example_uns"],
}
## VIASH END

################################################################################
# FUNCTIONS
################################################################################


def expand_slots(adata_base, adata_insert, slots):
    """
    Expand '__auto__' tokens in slots to automatically include items not in base.

    Parameters
    ----------
    adata_base : AnnData
        The base AnnData object
    adata_insert : AnnData
        The insert AnnData object
    slots : dict
        Dictionary with keys for each slot type containing lists of item names

    Returns
    -------
    dict
        Updated slots dictionary with '__auto__' tokens expanded
    """

    for slot in slots.keys():
        slot_items = set(slots[slot])
        if "__auto__" in slot_items:
            print(f"Expanding '__auto__' for slot '{slot}'...", flush=True)
            if slot == "obs":
                base_items = list(adata_base.obs.columns)
                insert_items = list(adata_insert.obs.columns)
            elif slot == "var":
                base_items = list(adata_base.var.columns)
                insert_items = list(adata_insert.var.columns)
            elif slot == "layers":
                base_items = list(adata_base.layers.keys())
                insert_items = list(adata_insert.layers.keys())
            elif slot == "obsm":
                base_items = list(adata_base.obsm.keys())
                insert_items = list(adata_insert.obsm.keys())
            elif slot == "varm":
                base_items = list(adata_base.varm.keys())
                insert_items = list(adata_insert.varm.keys())
            elif slot == "obsp":
                base_items = list(adata_base.obsp.keys())
                insert_items = list(adata_insert.obsp.keys())
            elif slot == "varp":
                base_items = list(adata_base.varp.keys())
                insert_items = list(adata_insert.varp.keys())
            elif slot == "uns":
                base_items = list(adata_base.uns.keys())
                insert_items = list(adata_insert.uns.keys())
            else:
                warn(
                    f"Unknown slot '{slot}', skipping __auto__ expansion.", stacklevel=2
                )
                continue

            new_items = [item for item in insert_items if item not in base_items]

            print(f"  Adding new items: {new_items}", flush=True)
            slots[slot] = list(slot_items - {"__auto__"}) + new_items

    return slots


def make_empty_series(source, length):
    """Create an empty Series with a nullable dtype compatible with source."""

    dtype = source.dtype

    if isinstance(dtype, pd.CategoricalDtype):
        return pd.Series(
            pd.Categorical(
                [pd.NA] * length,
                categories=source.cat.categories,
                ordered=source.cat.ordered,
            )
        )

    if pd.api.types.is_bool_dtype(dtype):
        return pd.Series(pd.array([pd.NA] * length, dtype="boolean"))

    if pd.api.types.is_integer_dtype(dtype):
        return pd.Series(pd.array([pd.NA] * length, dtype="Int64"))

    if pd.api.types.is_string_dtype(dtype) or dtype is object:
        return pd.Series(pd.array([pd.NA] * length, dtype="string"))

    return pd.Series(pd.array([np.nan] * length, dtype=dtype))


def insert_dataframe_column(
    df_base,
    df_insert,
    col,
    insert_indices_in_base,
    insert_indices_in_insert,
    df_name,
):
    """
    Insert a column from insert DataFrame into base DataFrame.

    Parameters
    ----------
    df_base : pd.DataFrame
        The base DataFrame (e.g., adata_base.obs)
    df_insert : pd.DataFrame
        The insert DataFrame (e.g., adata_insert.obs)
    col : str
        The column name to insert
    insert_indices_in_base : np.ndarray
        Indices in base where insert items map to
    insert_indices_in_insert : np.ndarray
        Indices in insert
    df_name : str
        Name of the DataFrame (obs/var) for printing
    """
    if col in df_base.columns:
        warn(
            f"  {df_name} column '{col}' already exists in base, content will be overwritten.",
            stacklevel=3,
        )

    if col in df_insert.columns:
        print(f"  Inserting {df_name} column: {col}", flush=True)

        col_data = df_insert[col]
        if pd.api.types.is_string_dtype(col_data.dtype):
            warn(
                f"  {df_name} column '{col}' has string dtype, converting to categorical",
                stacklevel=3,
            )
            col_data = col_data.astype("category")

        # Create empty series and insert data
        df_base[col] = make_empty_series(col_data, len(df_base))
        df_base.iloc[insert_indices_in_base, df_base.columns.get_loc(col)] = (
            col_data.iloc[insert_indices_in_insert].array
        )
    else:
        warn(
            f"  {df_name} column '{col}' not found in insert, skipping.",
            stacklevel=3,
        )


def insert_slot_item(
    slot_data_base,
    slot_data_insert,
    key,
    insert_indices_in_base,
    insert_indices_in_insert,
    index_names,
    slot_name,
    base_shape=None,
    missing_base_row_indices=None,
    missing_base_col_indices=None,
):
    """
    Generic function to insert data (dense, sparse, or DataFrame) from insert into base.

    Parameters
    ----------
    slot_data_base : dict-like
        The slot data container in base (e.g., adata_base.obsm)
    slot_data_insert : dict-like
        The slot data container in insert (e.g., adata_insert.obsm)
    key : str
        The key for the item to insert
    insert_indices_in_base : np.ndarray
        Indices in base where insert items map to
    insert_indices_in_insert : np.ndarray
        Indices in insert
    index_names : pd.Index
        Index labels for the output (obs_names or var_names)
    slot_name : str
        Name of the slot (e.g. obsm/varm/obsp/varp) for printing
    base_shape : tuple, optional
        Shape for the base container. If None, data is treated as dense.
        For symmetric matrices (obsp/varp), pass (n, n).
    missing_base_row_indices : np.ndarray, optional
        Indices in base rows not represented in insert.
    missing_base_col_indices : np.ndarray, optional
        Indices in base columns not represented in insert.
    """
    if key in slot_data_base and key in slot_data_insert:
        warn(
            f"  {slot_name} '{key}' already exists in base, content will be overwritten.",
            stacklevel=3,
        )

    if key in slot_data_insert:
        print(f"  Inserting {slot_name}: {key}", flush=True)
        data_insert = slot_data_insert[key]
        is_sparse = sp.issparse(data_insert)
        is_dataframe = isinstance(data_insert, pd.DataFrame)

        if base_shape is not None:
            # Handle square arrays (e.g. obsp/varp)
            insert_data = data_insert[
                np.ix_(insert_indices_in_insert, insert_indices_in_insert)
            ]
            if is_sparse:
                slot_data_base[key] = build_sparse_matrix_from_submatrix(
                    insert_data,
                    insert_indices_in_base,
                    insert_indices_in_base,
                    base_shape,
                    missing_base_row_indices=missing_base_row_indices,
                    missing_base_col_indices=missing_base_col_indices,
                    materialize_missing_nan=False,
                )
            else:
                base_array = np.full(base_shape, np.nan)
                base_array[np.ix_(insert_indices_in_base, insert_indices_in_base)] = (
                    np.asarray(insert_data)
                )
                slot_data_base[key] = base_array
        else:
            if is_dataframe:
                # Handle DataFrames column-by-column to preserve dtypes
                df_base = pd.DataFrame(index=index_names)
                for col in data_insert.columns:
                    insert_dataframe_column(
                        df_base,
                        data_insert,
                        col,
                        insert_indices_in_base,
                        insert_indices_in_insert,
                        slot_name,
                    )
                slot_data_base[key] = df_base
            elif is_sparse:
                insert_submatrix = data_insert[insert_indices_in_insert, :]
                n_cols = insert_submatrix.shape[1]
                base_col_indices = np.arange(n_cols)
                slot_data_base[key] = build_sparse_matrix_from_submatrix(
                    insert_submatrix,
                    insert_indices_in_base,
                    base_col_indices,
                    (len(index_names), n_cols),
                    missing_base_row_indices=missing_base_row_indices,
                    materialize_missing_nan=True,
                )
            else:
                data_array = np.asarray(data_insert)
                insert_submatrix = data_array[insert_indices_in_insert, :]
                base_array = np.full(
                    (len(index_names), insert_submatrix.shape[1]), np.nan
                )
                base_array[insert_indices_in_base, :] = insert_submatrix
                slot_data_base[key] = base_array
    else:
        warn(
            f"  {slot_name} '{key}' not found in insert, skipping.",
            stacklevel=3,
        )


def build_sparse_matrix_from_submatrix(
    data_submatrix,
    base_row_indices,
    base_col_indices,
    base_shape,
    missing_base_row_indices=None,
    missing_base_col_indices=None,
    materialize_missing_nan=True,
):
    """Build sparse output in base coordinates from a sparse insert submatrix.

    Parameters
    ----------
    data_submatrix : sparse matrix
        The sparse insert submatrix to reindex.
    base_row_indices : np.ndarray
        Row indices in base coordinates.
    base_col_indices : np.ndarray
        Column indices in base coordinates.
    base_shape : tuple
        Output shape.
    missing_base_row_indices : np.ndarray, optional
        Pre-computed indices of missing rows in base.
    missing_base_col_indices : np.ndarray, optional
        Pre-computed indices of missing columns in base.
    materialize_missing_nan : bool, default True
        If True, explicitly store NaN for unmapped base rows/columns.
        If False, missing regions are implicit (treated as 0).
    """

    coo = data_submatrix.tocoo()

    if coo.nnz > 0:
        mapped_rows = base_row_indices[coo.row]
        mapped_cols = base_col_indices[coo.col]
        data_values = coo.data
    else:
        mapped_rows = np.array([], dtype=np.intp)
        mapped_cols = np.array([], dtype=np.intp)
        data_values = np.array([], dtype=float)

    if not materialize_missing_nan:
        return sp.csr_matrix(
            (data_values, (mapped_rows, mapped_cols)),
            shape=base_shape,
        )

    n_rows, n_cols = base_shape

    if missing_base_row_indices is None:
        present_row_mask = np.zeros(n_rows, dtype=bool)
        present_row_mask[np.unique(base_row_indices)] = True
        missing_rows = np.flatnonzero(~present_row_mask)
    else:
        missing_rows = np.asarray(missing_base_row_indices, dtype=np.intp)
        present_row_mask = np.ones(n_rows, dtype=bool)
        present_row_mask[missing_rows] = False
    present_rows = np.flatnonzero(present_row_mask)

    if missing_base_col_indices is None:
        present_col_mask = np.zeros(n_cols, dtype=bool)
        present_col_mask[np.unique(base_col_indices)] = True
        missing_cols = np.flatnonzero(~present_col_mask)
    else:
        missing_cols = np.asarray(missing_base_col_indices, dtype=np.intp)

    nan_rows_a = (
        np.repeat(missing_rows, n_cols)
        if missing_rows.size
        else np.array([], dtype=np.intp)
    )
    nan_cols_a = (
        np.tile(np.arange(n_cols, dtype=np.intp), missing_rows.size)
        if missing_rows.size
        else np.array([], dtype=np.intp)
    )

    nan_rows_b = (
        np.repeat(present_rows, missing_cols.size)
        if missing_cols.size
        else np.array([], dtype=np.intp)
    )
    nan_cols_b = (
        np.tile(missing_cols, present_rows.size)
        if missing_cols.size
        else np.array([], dtype=np.intp)
    )

    nan_count = nan_rows_a.size + nan_rows_b.size
    if nan_count == 0:
        return sp.csr_matrix(
            (data_values, (mapped_rows, mapped_cols)), shape=base_shape
        )

    data_values = np.asarray(data_values, dtype=float)
    rows = np.concatenate([mapped_rows, nan_rows_a, nan_rows_b])
    cols = np.concatenate([mapped_cols, nan_cols_a, nan_cols_b])
    vals = np.concatenate([data_values, np.full(nan_count, np.nan)])

    return sp.csr_matrix((vals, (rows, cols)), shape=base_shape)


def insert_partial_items(
    mdata_base, mdata_insert, modality_base, modality_insert, slots
):
    """
    Insert items from mdata_insert into mdata_base.

    Parameters
    ----------
    mdata_base : MuData
        The base MuData object to insert into
    mdata_insert : MuData
        The MuData object to insert from
    modality_base : str
        The modality name in the base MuData
    modality_insert : str
        The modality name in the insert MuData
    slots : dict
        Dictionary with keys 'obs', 'var', 'layers', 'obsm', 'varm', 'obsp', 'varp', 'uns'
        Each key contains a list of slot names to insert
    """

    if modality_base not in mdata_base.mod:
        raise ValueError(f"Modality '{modality_base}' not found in base MuData")
    adata_base = mdata_base[modality_base]
    if modality_insert not in mdata_insert.mod:
        raise ValueError(f"Modality '{modality_insert}' not found in insert MuData")
    adata_insert = mdata_insert[modality_insert]

    # Check that all cells in insert are present in base
    if not adata_insert.obs_names.isin(adata_base.obs_names).all():
        missing_cells = adata_insert.obs_names[
            ~adata_insert.obs_names.isin(adata_base.obs_names)
        ]
        raise ValueError(
            f"{len(missing_cells)} cells in insert not found in base: {missing_cells}"
        )

    # Check that all variables in insert are present in base
    if not adata_insert.var_names.isin(adata_base.var_names).all():
        missing_vars = adata_insert.var_names[
            ~adata_insert.var_names.isin(adata_base.var_names)
        ]
        raise ValueError(
            f"{len(missing_vars)} variables in insert not found in base: {missing_vars}"
        )

    # Get indices of cells in insert object within the base object
    common_obs_names = adata_base.obs_names.intersection(adata_insert.obs_names)
    insert_indices_in_base = adata_base.obs_names.get_indexer(common_obs_names)
    insert_indices_in_insert = adata_insert.obs_names.get_indexer(common_obs_names)
    present_obs_mask = np.zeros(adata_base.n_obs, dtype=bool)
    present_obs_mask[insert_indices_in_base] = True
    missing_obs_indices_in_base = np.flatnonzero(~present_obs_mask)

    # Get indices of variables in insert object within the base object
    common_var_names = adata_base.var_names.intersection(adata_insert.var_names)
    insert_var_indices_in_base = adata_base.var_names.get_indexer(common_var_names)
    insert_var_indices_in_insert = adata_insert.var_names.get_indexer(common_var_names)
    present_var_mask = np.zeros(adata_base.n_vars, dtype=bool)
    present_var_mask[insert_var_indices_in_base] = True
    missing_var_indices_in_base = np.flatnonzero(~present_var_mask)

    print("Inserting obs columns...", flush=True)
    for col in slots["obs"]:
        insert_dataframe_column(
            adata_base.obs,
            adata_insert.obs,
            col,
            insert_indices_in_base,
            insert_indices_in_insert,
            "obs",
        )

    print("Inserting var columns...", flush=True)
    for col in slots["var"]:
        insert_dataframe_column(
            adata_base.var,
            adata_insert.var,
            col,
            insert_var_indices_in_base,
            insert_var_indices_in_insert,
            "var",
        )

    print("Inserting layers...", flush=True)
    for layer in slots["layers"]:
        if layer in adata_base.layers and layer in adata_insert.layers:
            warn(
                f"  Layer '{layer}' already exists in base, content will be overwritten.",
                stacklevel=2,
            )
        if layer in adata_insert.layers:
            print(f"  Inserting layer: {layer}", flush=True)
            insert_layer_data = adata_insert.layers[layer][
                np.ix_(insert_indices_in_insert, insert_var_indices_in_insert)
            ]
            if sp.issparse(insert_layer_data):
                adata_base.layers[layer] = build_sparse_matrix_from_submatrix(
                    insert_layer_data,
                    insert_indices_in_base,
                    insert_var_indices_in_base,
                    (adata_base.n_obs, adata_base.n_vars),
                    missing_base_row_indices=missing_obs_indices_in_base,
                    missing_base_col_indices=missing_var_indices_in_base,
                    materialize_missing_nan=True,
                )
            else:
                adata_base.layers[layer] = np.full(
                    (adata_base.n_obs, adata_base.n_vars), np.nan
                )
                adata_base.layers[layer][
                    np.ix_(insert_indices_in_base, insert_var_indices_in_base)
                ] = np.asarray(insert_layer_data)
        else:
            warn(f"  Layer '{layer}' not found in insert, skipping.", stacklevel=2)

    print("Inserting obsm...", flush=True)
    for key in slots["obsm"]:
        insert_slot_item(
            adata_base.obsm,
            adata_insert.obsm,
            key,
            insert_indices_in_base,
            insert_indices_in_insert,
            adata_base.obs_names,
            "obsm",
            missing_base_row_indices=missing_obs_indices_in_base,
        )

    print("Inserting varm...", flush=True)
    for key in slots["varm"]:
        insert_slot_item(
            adata_base.varm,
            adata_insert.varm,
            key,
            insert_var_indices_in_base,
            insert_var_indices_in_insert,
            adata_base.var_names,
            "varm",
            missing_base_row_indices=missing_var_indices_in_base,
        )

    print("Inserting obsp...", flush=True)
    for key in slots["obsp"]:
        insert_slot_item(
            adata_base.obsp,
            adata_insert.obsp,
            key,
            insert_indices_in_base,
            insert_indices_in_insert,
            adata_base.obs_names,
            "obsp",
            base_shape=(adata_base.n_obs, adata_base.n_obs),
            missing_base_row_indices=missing_obs_indices_in_base,
            missing_base_col_indices=missing_obs_indices_in_base,
        )

    print("Inserting varp...", flush=True)
    for key in slots["varp"]:
        insert_slot_item(
            adata_base.varp,
            adata_insert.varp,
            key,
            insert_var_indices_in_base,
            insert_var_indices_in_insert,
            adata_base.var_names,
            "varp",
            base_shape=(adata_base.n_vars, adata_base.n_vars),
            missing_base_row_indices=missing_var_indices_in_base,
            missing_base_col_indices=missing_var_indices_in_base,
        )

    print("Inserting uns...", flush=True)
    for key in slots["uns"]:
        if key in adata_base.uns and key in adata_insert.uns:
            warn(
                f"  Uns '{key}' already exists in base, content will be overwritten.",
                stacklevel=2,
            )
        if key in adata_insert.uns:
            print(f"  Inserting uns: {key}", flush=True)
            adata_base.uns[key] = adata_insert.uns[key]
        else:
            warn(f"  Uns '{key}' not found in insert, skipping.", stacklevel=2)

    return mdata_base


################################################################################
# MAIN
################################################################################


def main(par):
    print(
        f"====== Insert partial items (mudata v{mudata.__version__}) ======", flush=True
    )

    print(f"\n>>> Loading base MuData from '{par['input_base']}'...", flush=True)
    mdata_base = mudata.read_h5mu(par["input_base"])
    print(mdata_base, flush=True)
    print(f"Modality '{par['modality_base']}':", flush=True)

    print(f"\n>>> Loading insert MuData from '{par['input_insert']}'...", flush=True)
    mdata_insert = mudata.read_h5mu(par["input_insert"])
    print(mdata_insert, flush=True)
    print(f"Modality '{par['modality_insert']}':", flush=True)

    print("\n>>> Expanding slots...", flush=True)
    slots = {
        "obs": par["obs"] or [],
        "var": par["var"] or [],
        "layers": par["layers"] or [],
        "obsm": par["obsm"] or [],
        "varm": par["varm"] or [],
        "obsp": par["obsp"] or [],
        "varp": par["varp"] or [],
        "uns": par["uns"] or [],
    }
    print(f"Initial slots: {slots}", flush=True)
    slots = expand_slots(
        mdata_base[par["modality_base"]], mdata_insert[par["modality_insert"]], slots
    )
    print(f"Expanded slots: {slots}", flush=True)

    print("\n>>> Inserting items...", flush=True)
    mdata_base = insert_partial_items(
        mdata_base, mdata_insert, par["modality_base"], par["modality_insert"], slots
    )

    print(f"\n>>> Writing output to '{par['output']}'...", flush=True)
    print(mdata_base, flush=True)
    print(f"Modality '{par['modality_base']}':", flush=True)
    print(mdata_base[par["modality_base"]], flush=True)
    mdata_base.write_h5mu(par["output"])

    print("\n>>> Done!\n")


if __name__ == "__main__":
    sys.exit(main(par))
