import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp
from mudata import read_h5ad

## VIASH START
par = {
    "input": "input.h5mu",
    "modality": "rna",
    "source": "source.h5mu",
    "source_modality": None,
    "source_var_gene_names": None,
    "obs": None,
    "var": None,
    "layers": None,
    "obsm": None,
    "varm": None,
    "obsp": None,
    "varp": None,
    "uns": None,
    "allow_overwrite": False,
    "allow_partial": False,
    "output": "output.h5mu",
    "output_compression": None,
}
meta = {"resources_dir": "src/utils/"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger  # noqa: E402
from compress_h5mu import write_h5ad_to_h5mu_with_compression  # noqa: E402

logger = setup_logger()

# Slots grouped by which axes they index. layers is special (both axes).
_OBS_ONLY_SLOTS = ("obs", "obsm")
_VAR_ONLY_SLOTS = ("var", "varm")
_OBS_SQUARE_SLOTS = ("obsp",)
_VAR_SQUARE_SLOTS = ("varp",)
_LAYERS_SLOT = "layers"
_UNS_SLOT = "uns"
_ALL_SLOTS = (
    "obs",
    "var",
    "layers",
    "obsm",
    "varm",
    "obsp",
    "varp",
    "uns",
)
# Which slots index each axis. layers indexes both, so it appears under both.
_AXIS_SLOTS = {
    "obs": (*_OBS_ONLY_SLOTS, *_OBS_SQUARE_SLOTS, _LAYERS_SLOT),
    "var": (*_VAR_ONLY_SLOTS, *_VAR_SQUARE_SLOTS, _LAYERS_SLOT),
}
# Nullable dtype to hold an unmatched position, per dtype of the source column.
# Checked in order; anything unmatched keeps its own dtype and is NaN-filled.
_EMPTY_DTYPES = (
    (pd.api.types.is_bool_dtype, "boolean"),
    (pd.api.types.is_integer_dtype, "Int64"),
    (pd.api.types.is_string_dtype, "string"),
)


def slot_items(mod, slot):
    """List the item names held in ``slot`` of modality ``mod``."""
    value = getattr(mod, slot)
    return list(value.columns if slot in ("obs", "var") else value.keys())


def resolve_slots(source_mod, input_mod, slots):
    """Resolve each slot to the items to copy. A slot left unspecified copies
    every item present in the source but missing from the input. Mutates and
    returns ``slots`` in place."""
    for slot, keys in slots.items():
        if keys:
            continue
        input_items = slot_items(input_mod, slot)
        slots[slot] = [
            item for item in slot_items(source_mod, slot) if item not in input_items
        ]
        if slots[slot]:
            logger.info("No .%s specified; copying %s", slot, slots[slot])
    return slots


def make_empty_series(source, length):
    """Create an empty Series with a nullable dtype compatible with ``source``."""
    dtype = source.dtype
    if isinstance(dtype, pd.CategoricalDtype):
        return pd.Series(
            pd.Categorical(
                [pd.NA] * length,
                categories=source.cat.categories,
                ordered=source.cat.ordered,
            )
        )
    for matches, empty_dtype in _EMPTY_DTYPES:
        if matches(dtype):
            return pd.Series(pd.array([pd.NA] * length, dtype=empty_dtype))
    return pd.Series(pd.array([np.nan] * length, dtype=dtype))


def reindex_series(source_series, source_to_target_idx, n_target):
    """Build a target-length Series with ``source_series`` placed at the given
    target positions; remaining positions hold the dtype's missing value."""
    source_col = source_series
    if pd.api.types.is_string_dtype(source_col.dtype):
        # Strings get categorical-ised so missing values can be represented
        # naturally; aligns with the convention used elsewhere in the codebase.
        source_col = source_col.astype("category")
    out = make_empty_series(source_col, n_target)
    out.iloc[source_to_target_idx] = source_col.array
    return out


def reindex_dense_axis(data, source_to_target_idx, n_target, axis):
    """Reindex a dense array along one axis, NaN-filling unmapped positions."""
    data = np.asarray(data)
    if axis == 0:
        out = np.full((n_target,) + data.shape[1:], np.nan, dtype=float)
        out[source_to_target_idx] = data
    else:
        out = np.full(
            data.shape[:1] + (n_target,) + data.shape[2:], np.nan, dtype=float
        )
        out[:, source_to_target_idx] = data
    return out


def reindex_sparse_rows(data, source_to_target_rows, n_target_rows):
    coo = data.tocoo()
    mapped_rows = source_to_target_rows[coo.row]
    return sp.csr_matrix(
        (coo.data, (mapped_rows, coo.col)),
        shape=(n_target_rows, coo.shape[1]),
    )


def reindex_sparse_cols(data, source_to_target_cols, n_target_cols):
    coo = data.tocoo()
    mapped_cols = source_to_target_cols[coo.col]
    return sp.csr_matrix(
        (coo.data, (coo.row, mapped_cols)),
        shape=(coo.shape[0], n_target_cols),
    )


def reindex_2d(data, src_to_tgt_rows, src_to_tgt_cols, n_rows, n_cols):
    """Reindex both axes of an array (or square sparse matrix)."""
    if sp.issparse(data):
        coo = data.tocoo()
        return sp.csr_matrix(
            (coo.data, (src_to_tgt_rows[coo.row], src_to_tgt_cols[coo.col])),
            shape=(n_rows, n_cols),
        )
    data = np.asarray(data)
    out = np.full((n_rows, n_cols), np.nan, dtype=float)
    out[np.ix_(src_to_tgt_rows, src_to_tgt_cols)] = data
    return out


def reindex_obsm_value(value, src_to_tgt_obs, target_axis_index):
    """Reindex a single .obsm/.varm entry along its first axis."""
    n_target = len(target_axis_index)
    if isinstance(value, pd.DataFrame):
        out = pd.DataFrame(index=target_axis_index.copy())
        for col in value.columns:
            new_col = reindex_series(value[col], src_to_tgt_obs, n_target)
            out[col] = new_col.values
        return out
    if sp.issparse(value):
        return reindex_sparse_rows(value, src_to_tgt_obs, n_target)
    return reindex_dense_axis(value, src_to_tgt_obs, n_target, axis=0)


def validate_axis_containment(source_idx, target_idx, axis, allow_partial):
    """Return source-to-target position mapping; raise if the source has names not
    in the target, or (unless ``allow_partial``) the sets differ."""
    src_to_tgt = target_idx.get_indexer(source_idx)
    missing_in_target = src_to_tgt < 0
    if missing_in_target.any():
        missing = source_idx[missing_in_target].tolist()
        raise ValueError(
            f"Source .{axis}_names contains {missing_in_target.sum()} entries "
            f"not present in input .{axis}_names (first 5: {missing[:5]}). "
            f"Source must be a subset of input."
        )
    if not allow_partial and len(source_idx) != len(target_idx):
        raise ValueError(
            f"Source and input .{axis}_names differ in size "
            f"(source: {len(source_idx)}, input: {len(target_idx)}). "
            f"Use --allow_partial to permit a strict subset and NaN-fill the "
            f"remaining input rows/columns."
        )
    return src_to_tgt


def copy_dataframe_columns(source_df, target_df, columns, src_to_tgt, axis_name):
    """Copy named columns from ``source_df`` into ``target_df``, reindexing
    rows from source positions to target positions."""
    n_target = len(target_df)
    for col in columns:
        new_series = reindex_series(source_df[col], src_to_tgt, n_target)
        # Use positional assignment to avoid pandas trying to align on labels.
        target_df[col] = new_series.values
    logger.info("Copied .%s columns: %s", axis_name, columns)


def main():
    source_modality = par["source_modality"] or par["modality"]

    logger.info(
        "Reading modality '%s' from source file '%s'",
        source_modality,
        par["source"],
    )
    try:
        source_mod = read_h5ad(par["source"], mod=source_modality)
    except KeyError:
        raise ValueError(
            f"Modality '{source_modality}' does not exist in source "
            f"file '{par['source']}'."
        )

    logger.info(
        "Reading modality '%s' from input file '%s'",
        par["modality"],
        par["input"],
    )
    try:
        input_mod = read_h5ad(par["input"], mod=par["modality"])
    except KeyError:
        raise ValueError(
            f"Modality '{par['modality']}' does not exist in input file "
            f"'{par['input']}'."
        )

    if par["source_var_gene_names"]:
        col = par["source_var_gene_names"]
        if col not in source_mod.var.columns:
            raise ValueError(
                f"--source_var_gene_names '{col}' not found in source .var "
                f"columns: {list(source_mod.var.columns)}"
            )
        logger.info("Overriding source .var index with column '%s'", col)
        source_mod.var_names = pd.Index(source_mod.var[col].astype(str))

    slots = {name: par[name] or [] for name in _ALL_SLOTS}
    slots = resolve_slots(source_mod, input_mod, slots)
    requested = {k: v for k, v in slots.items() if v}
    if not requested:
        logger.info("Nothing to copy; writing input unchanged.")

    src_to_tgt = dict.fromkeys(_AXIS_SLOTS)
    for axis, axis_slots in _AXIS_SLOTS.items():
        if not any(slot in requested for slot in axis_slots):
            continue
        src_to_tgt[axis] = validate_axis_containment(
            getattr(source_mod, f"{axis}_names"),
            getattr(input_mod, f"{axis}_names"),
            axis,
            par["allow_partial"],
        )

    # Check existence in source and overwrite policy in input for every slot.
    for slot, keys in requested.items():
        available_source = set(slot_items(source_mod, slot))
        available_input = set(slot_items(input_mod, slot))
        missing = [k for k in keys if k not in available_source]
        if missing:
            raise ValueError(
                f"The following .{slot} keys were not found in source modality "
                f"'{source_modality}': {missing}"
            )
        existing = [k for k in keys if k in available_input]
        if existing and not par["allow_overwrite"]:
            raise ValueError(
                f"The following .{slot} keys already exist in input modality "
                f"'{par['modality']}': {existing}. Use --allow_overwrite to "
                f"overwrite them."
            )
        if existing:
            logger.warning("Overwriting existing .%s keys: %s", slot, existing)

    n_input_obs = input_mod.n_obs
    n_input_var = input_mod.n_vars

    # .obs / .var
    for slot in ("obs", "var"):
        if slot in requested:
            copy_dataframe_columns(
                getattr(source_mod, slot),
                getattr(input_mod, slot),
                requested[slot],
                src_to_tgt[slot],
                slot,
            )

    # .obsm / .varm
    for slot, axis, input_axis_index in (
        ("obsm", "obs", input_mod.obs_names),
        ("varm", "var", input_mod.var_names),
    ):
        for key in requested.get(slot, []):
            value = getattr(source_mod, slot)[key]
            getattr(input_mod, slot)[key] = reindex_obsm_value(
                value, src_to_tgt[axis], input_axis_index
            )
            logger.info("Copied .%s key: %s", slot, key)

    # .obsp / .varp (square)
    for slot, axis, n_input in (
        ("obsp", "obs", n_input_obs),
        ("varp", "var", n_input_var),
    ):
        for key in requested.get(slot, []):
            value = getattr(source_mod, slot)[key]
            getattr(input_mod, slot)[key] = reindex_2d(
                value, src_to_tgt[axis], src_to_tgt[axis], n_input, n_input
            )
            logger.info("Copied .%s key: %s", slot, key)

    # .layers (both axes)
    for key in requested.get("layers", []):
        value = source_mod.layers[key]
        input_mod.layers[key] = reindex_2d(
            value, src_to_tgt["obs"], src_to_tgt["var"], n_input_obs, n_input_var
        )
        logger.info("Copied .layers key: %s", key)

    # .uns (no reindexing)
    for key in requested.get("uns", []):
        input_mod.uns[key] = source_mod.uns[key]
        logger.info("Copied .uns key: %s", key)

    logger.info(
        "Writing output to '%s' with compression '%s'",
        par["output"],
        par["output_compression"],
    )
    write_h5ad_to_h5mu_with_compression(
        output_file=par["output"],
        h5mu=par["input"],
        modality_name=par["modality"],
        modality_data=input_mod,
        output_compression=par["output_compression"],
    )


if __name__ == "__main__":
    sys.exit(main())
