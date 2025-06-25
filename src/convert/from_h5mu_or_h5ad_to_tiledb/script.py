import sys
import anndata
import numpy as np
import tiledbsoma

# This import seems redundant, but it is required!
import tiledbsoma.io
import pandas as pd
from functools import singledispatch
from collections.abc import Mapping
import h5py
from shutil import copy2, rmtree
from collections.abc import Callable
from functools import partial
from pathlib import Path
from tempfile import NamedTemporaryFile
import pyarrow as pa
from collections import defaultdict
import json

anndata.settings.allow_write_nullable_strings = True


## VIASH START
par = {
    "input": "/home/di/code/openpipelines-multisample/pmbc_process_samples.h5mu",
    # modality names
    "rna_modality": "rna",
    "rna_modality_output": "rna",
    "prot_modality": "prot",
    "prot_modality_output": "prot",
    # rna var
    "rna_var_gene_names_input": "gene_symbol",
    "rna_var_gene_names_output": "gene_foo",
    "rna_var_index_name_output": "rna_index",
    # rna layers
    "rna_raw_layer_input": "X",
    "rna_raw_layer_output": "X",
    "rna_normalized_layer_input": "log_normalized",
    "rna_normalized_layer_output": "log_normalized",
    # prot layers
    "prot_raw_layer_input": "X",
    "prot_raw_layer_output": "X",
    "prot_normalized_layer_input": "clr",
    "prot_normalized_layer_output": "log_normalized",
    # prot var
    "prot_var_index_name_output": "prot_index",
    # observation index name
    "obs_index_name_output": "cell_id",
    # TileDB output directory
    "tiledb_dir": "output",
}


meta = {
    "name": "from_h5mu_or_h5ad_to_tiledb",
    "resources_dir": "src/utils/",
    "temp_dir": "/tmp",
}
## VIASH END
sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

LAYER_ARGUMENTS = (
    "rna_raw_layer_input",
    "rna_raw_layer_output",
    "rna_normalized_layer_input",
    "rna_normalized_layer_output",
    "prot_raw_layer_input",
    "prot_raw_layer_output",
    "prot_normalized_layer_input",
    "prot_normalized_layer_output",
)


multidim_column_indices = {}


def _convert_to_array(mod_name, _, group):
    """
    tiledb SOMA does not support pandas dataframes in varm and obsm
    Convert those to structured arrays and store column index as metadata
    """
    entries = group.keys()
    for entry in entries:
        encoding = group[entry].attrs["encoding-type"]
        if encoding in ("array", "csr_matrix", "csc_matrix"):
            return
        if encoding in {
            "dataframe",
            "nullable-integer",
            "nullable-boolean",
            "nullable-string-array",
        }:
            df = anndata.io.read_elem(group[entry])
            data = df.to_numpy()
            anndata.io.write_elem(group, entry, data)
            index_to_write = df.columns
            # No need to store write the index to tileDB in the case its just the integer index
            if not isinstance(index_to_write, pd.RangeIndex):
                multidim_column_indices.setdefault(mod_name, {}).setdefault(
                    group.name.strip("/"), {}
                ).update({entry: index_to_write.to_list()})
            return
        logger.info(
            "Deleting %s because it is of encoding %s, which is not supported.",
            f"{group.name}/{entry}",
            encoding,
        )
        del group[entry]


def _read_obs(anndata_file):
    with h5py.File(anndata_file, "r") as open_rna:
        return _to_pandas_nullable(anndata.io.read_elem(open_rna["obs"]))


def _write_obs(anndata_file, obs_df):
    with h5py.File(anndata_file, "r+") as open_h5:
        anndata.io.write_elem(open_h5, "obs", obs_df)


def _log_arguments(function_obj, arg_dict):
    """
    Format a dictionairy of arguments into a string that is put into the script logs.
    """
    args_str = [f"\t{param}: {param_val}\n" for param, param_val in arg_dict.items()]
    logger.info(
        "Calling %s with arguments:\n%s",
        function_obj.__name__,
        "".join(args_str).rstrip(),
    )


def _to_pandas_nullable(table_or_df: pa.Table | pd.DataFrame) -> pd.DataFrame:
    """
    Convert pyarrow Table or pandas dataframe to pandas dataframe that uses nullable data types.

    This function takes the following into account:
        1. AnnData supports writing BooleanDtype, IntegerDtype and StringDtype (which uses pd.NA)
        2. AnnData does *not* support FloatingDtype
        3. AnnData does support numpy-stype floating and integer dtypes (which use np.nan).
        4. TileDB SOMA's `update_obs` does not allow casting floats to integers (potential information loss).
        5. When pandas's `convert_dtypes` convert_integer argument is set to True,
           floats will be converted to integer if they can be 'faithfully' casted.

    Giving the previous constraints, it is *not* possible to only convert integers using `convert_types` because
    it might cause some float columns to be converted to integers (which is incompatible with tiledbsoma).
    Therefore, this function will leave integers and floats as-is.

    Based on the conversion provided by this function, a suitable na representation (pd.NA or np.NA) for a column can
    be chosen using `pd.api.types.is_extension_array_dtype(column_dtype)`
    """
    try:
        df = table_or_df.to_pandas(
            types_mapper=defaultdict(None, {pa.large_string(): pd.StringDtype()}).get
        )
    except AttributeError:
        df = table_or_df.convert_dtypes(
            infer_objects=False,
            convert_string=True,
            convert_integer=False,
            convert_boolean=False,
            convert_floating=False,
        )
    return df.convert_dtypes(
        infer_objects=False,
        convert_string=False,
        convert_integer=False,
        convert_boolean=True,
        convert_floating=False,
    )


def _get_temp_h5ad():
    return Path(NamedTemporaryFile(suffix=".h5ad", dir=meta["temp_dir"]).name)


def _h5mu_to_h5ad(h5mu_path, modality_name):
    """
    Create an anndata file from a modality in a mudata file by
    copying over the relevant h5 objects.
    """
    result = _get_temp_h5ad()
    with h5py.File(h5mu_path, "r") as open_anndata:
        with h5py.File(result, "w") as file_to_write:
            for h5_key in open_anndata["mod"][modality_name].keys():
                open_anndata.copy(
                    source=f"/mod/{modality_name}/{h5_key}", dest=file_to_write
                )
    return result


def _handle_layer_location(layer_name):
    """
    Handle the special case where 'X' is specified as layer
    this layer is not present in `<mod_name>/layers/X` but in `<mod_name>/X`.
    """
    prefix = "layers/" if layer_name != "X" else ""
    return f"{prefix}{layer_name}"


def _check_input_args(par, input, is_mudata):
    """
    Check the input arguments:
        * Arguments for protein modalities are not allowed when the input is not a MuData file.
        * Arguments for input layers may not point to the same layer
        * Output locations for layers may not be the same
        * Input modalities must not be the same
        * Input modalities must exist
    """
    if not is_mudata:
        for par_name in (
            "prot_modality",
            "prot_raw_layer_input",
            "prot_normalized_layer_input",
        ):
            if par[par_name]:
                raise ValueError(
                    f"'{par_name}' can not be used when the input is not a MuData file."
                )

    if par["prot_modality"]:
        for prot_arg in ("prot_raw_layer_input", "prot_raw_layer_output"):
            if not par[prot_arg]:
                ValueError(
                    f"When providing 'prot_modality', '{prot_arg}' must also be set."
                )

    NOT_SAME_VALUE = (
        ("rna_raw_layer_input", "rna_normalized_layer_input"),
        ("prot_raw_layer_input", "prot_normalized_layer_input"),
        ("rna_raw_layer_output", "rna_normalized_layer_output"),
        ("prot_raw_layer_output", "prot_normalized_layer_output"),
        ("rna_modality", "prot_modality"),
    )
    for layer1_arg, layer2_arg in NOT_SAME_VALUE:
        arg1_val, arg2_val = par[layer1_arg], par[layer2_arg]
        if (arg1_val == arg2_val) and arg1_val is not None:
            raise ValueError(
                f"The value for argument '{layer1_arg}' ({arg1_val}) must not be the same as for argument '{layer2_arg}' ({arg2_val})"
            )

    with h5py.File(input, "r") as open_h5:
        if "mod" in open_h5:
            available_modalities = tuple(open_h5["mod"].keys())
            for mod in (par["rna_modality"], par["prot_modality"]):
                if mod is not None and f"/mod/{mod}" not in open_h5:
                    raise ValueError(
                        f"Modality '{mod}' was not found in object, available modalities: {available_modalities}"
                    )


def _copy_layer(source: str, dest: str, input_h5ad: Path, mod_obj: h5py.Group):
    """
    Copy a h5py object from one location (source) to a new location (dest).
    The source location must be located in the input AnnData file, the destination
    will be written to a h5py Group object.
    """
    logger.info("Copying layer '%s' from '%s' to '%s'", source, input_h5ad, dest)
    if source == f"{mod_obj.name}/{dest}":
        logger.info("Layer is already at the correct location. Nothing to do.")
        return
    if not input_h5ad.is_file():
        raise FileNotFoundError(f"Could not link to {input_h5ad}: file does not exist.")
    with h5py.File(input_h5ad) as open_input:
        if source not in open_input:
            raise ValueError(f"Layer {source} was not found in {input_h5ad}.")
    if dest in mod_obj:
        logger.info(f"{dest} is already present. Removing before creating link.")
        del mod_obj[dest]
    # ExternalLink could have been used instead of 'copy', but it does not seem to work with tileDB
    with h5py.File(str(input_h5ad), "r") as open_input:
        open_input.copy(
            source=open_input[source], dest=mod_obj, name=dest, expand_external=True
        )
    logger.info("Copying complete.")


def _copy_column(src, dest, _, df_h5_obj):
    """
    Duplicate a column in a dataframe and give a new name.
    The input dataframe must be h5py object that is encoded by AnnData to store
    a pandas dataframe (e.g. .obs and .var).
    """
    logger.info("Copying column '%s' to '%s' for '%s'.", src, dest, df_h5_obj.name)
    columns = df_h5_obj.attrs["column-order"]
    if dest in columns:
        raise ValueError(f"Column {dest} already exists in {df_h5_obj.name}.")
    if src == dest:
        logger.info("Source and destination for the column are the same, not copying.")
        return None
    df_h5_obj.attrs["column-order"] = np.append(columns, dest)
    if src is None:
        # Use the index as source.
        src = df_h5_obj.attrs["_index"]
        # Both the index and columns are written as keys to df_h5_obj
        # An index name may be allowed to also be column name _only_ when
        # the column and the index have the same contents (this is an anndata requirement).
        # Here, this is always the case.
        if src in df_h5_obj.attrs["column-order"]:
            return None
    df_h5_obj.copy(src, dest, expand_refs=True, expand_soft=True)


def _set_index_name(name, _, df_h5_obj):
    """
    Set the index for a dataframe to an existing column.
    The input dataframe must be h5py object that is encoded by AnnData to store
    a pandas dataframe (e.g. .obs and .var).
    """
    logger.info("Setting index name of '%s' to '%s'.", df_h5_obj.name, name)
    original_index_name = df_h5_obj.attrs["_index"]
    # An index and a column may share a key in df_h5_obj
    # In this case we must not remove the original key from the object
    operator = (
        df_h5_obj.move
        if original_index_name not in df_h5_obj.attrs["column-order"]
        else df_h5_obj.copy
    )
    operator(original_index_name, name)
    df_h5_obj.attrs["_index"] = name
    logger.info("Done setting index name")


def _delete_all_keys(_, group_obj, exception=None):
    """
    Delete all keys from a h5py.Group object; which optional exceptions.
    """
    if exception is None:
        exception = ()
    logger.info(
        "Deleting all keys from '%s'%s.",
        group_obj.name,
        f" except {', '.join(exception)}" if exception else "",
    )
    keys_to_delete = set(group_obj.keys()) - set(exception)
    if not keys_to_delete:
        logger.info("No keys to delete.")
        return
    for key in keys_to_delete:
        del group_obj[key]
    logger.info("Done deleting keys %s.", ", ".join(keys_to_delete))


def _set_index_to_string_dtype(_, df_obj):
    """
    Change the index dataype for a dataframe to string when it is encoded as categorical.
    The dataframe must be provided as a h5py Group object that has been encoded by AnnData.
    """
    logger.info(
        "Checking if index of object '%s' is encoded as a categorical.", df_obj.name
    )
    index_col_name = df_obj.attrs["_index"]
    index_col = df_obj[index_col_name]
    logger.info("Column '%s' is being used as the index.", index_col_name)
    index_col_dtype = index_col.attrs["encoding-type"]
    if index_col_dtype == "categorical":
        logger.info(
            "Column '%s' is encoded as categorical, converting to string.",
            index_col_name,
        )
        is_ordered = bool(index_col.attrs.get("ordered", False))
        codes, categories = index_col["codes"], index_col["categories"]
        column = pd.Categorical.from_codes(codes, categories, is_ordered)
        column_str = column.astype(str).astype(object)
        del df_obj[index_col_name]
        df_obj.create_dataset(index_col_name, data=column_str)
        index_col.attrs["encoding-type"] = "string-array"
        logger.info("Conversion to string dtype complete.")
        return
    logger.info(
        "Column '%s' not encoded as categorical. Leaving as is.", index_col_name
    )


def _get_rna_conversion_specification(par):
    rna_spec = [
        partial(_copy_layer, par["rna_raw_layer_input"], par["rna_raw_layer_output"]),
        partial(
            _copy_layer,
            par["rna_normalized_layer_input"],
            par["rna_normalized_layer_output"],
        ),
        {
            "varm": partial(_convert_to_array, "rna"),
            "obsm": partial(_convert_to_array, "rna"),
            "uns": _delete_all_keys,
            "obsp": _delete_all_keys,
            "var": [
                partial(
                    _copy_column,
                    par["rna_var_gene_names_input"],
                    par["rna_var_gene_names_output"],
                ),
                partial(_set_index_name, par["rna_var_index_name_output"]),
                _set_index_to_string_dtype,
            ],
            "obs": partial(_set_index_name, par["obs_index_name_output"]),
            "layers": partial(
                _delete_all_keys,
                exception=(
                    par["rna_raw_layer_output"].removeprefix("layers/"),
                    par["rna_normalized_layer_output"].removeprefix("layers/"),
                ),
            ),
        },
    ]
    return rna_spec


def _get_prot_conversion_specification(par):
    prot_spec = [
        partial(_copy_layer, par["prot_raw_layer_input"], par["prot_raw_layer_output"]),
        partial(
            _copy_layer,
            par["prot_normalized_layer_input"],
            par["prot_normalized_layer_output"],
        ),
        {
            "varm": partial(_convert_to_array, "prot"),
            "obsm": partial(_convert_to_array, "prot"),
            "var": [
                partial(_set_index_name, par["prot_var_index_name_output"]),
                _set_index_to_string_dtype,
            ],
            "uns": _delete_all_keys,
            "obsp": _delete_all_keys,
            # "varm": _delete_all_keys,
            "layers": partial(
                _delete_all_keys,
                exception=(
                    par["prot_raw_layer_output"].removeprefix("layers/"),
                    par["prot_normalized_layer_output"].removeprefix("layers/"),
                ),
            ),
            "obs": partial(_set_index_name, par["obs_index_name_output"]),
        },
    ]

    return prot_spec


def _copy_anndata_and_convert_inplace(anndata_path, conversion_specification):
    converted_anndata = _get_temp_h5ad()
    copy2(anndata_path, converted_anndata)
    with h5py.File(converted_anndata, "r+") as open_prot_h5:
        apply_conversion(conversion_specification, anndata_path, open_prot_h5)
    return converted_anndata


@singledispatch
def apply_conversion(conversion_specification, input_h5ad, open_h5):
    raise NotImplementedError(
        f"Error: function 'apply_conversion' is not implemented for type '{type(conversion_specification)}'"
    )


@apply_conversion.register
def _(conversion_specification: list, input_h5ad, open_h5):
    for item in conversion_specification:
        apply_conversion(item, input_h5ad, open_h5)


@apply_conversion.register
def _(conversion_specifications: Mapping, input_h5ad, open_h5):
    for h5_key, conversion_spec in conversion_specifications.items():
        if hasattr(open_h5, h5_key):
            h5_obj = getattr(open_h5, h5_key)
        else:
            h5_obj = open_h5[h5_key]
        apply_conversion(conversion_spec, input_h5ad, h5_obj)


@apply_conversion.register
def _(conversion_specifications: Callable, input_h5ad, open_h5):
    conversion_specifications(input_h5ad, open_h5)


def _write_varm_obsm_indices(output_dir, mod_name, multidim_column_indices):
    # TileDB SOMA does not support column indices for .varm and obsm matrices
    # As a workaround, we will store these in the metadata slot
    for multidim_key in ("varm", "obsm"):
        collection_name = f"{str(output_dir)}/ms/{mod_name}/{multidim_key}"
        if tiledbsoma.Collection.exists(collection_name):
            with tiledbsoma.Collection.open(
                collection_name, "w"
            ) as multidim_collection:
                for item in multidim_collection:
                    index_to_write = (
                        multidim_column_indices.get(mod_name, {})
                        .get(multidim_key, {})
                        .get(item, None)
                    )
                    if index_to_write:
                        multidim_collection[item].metadata["column_index"] = json.dumps(
                            index_to_write
                        )


def _remove_directory(dir_path):
    """
    Check if a directory exists and remove it does.
    """
    if dir_path.exists():
        try:
            # Delete the directory and all its contents
            for item in dir_path.iterdir():
                if item.is_file() or item.is_symlink():
                    item.unlink()  # Removes files or symbolic links
                elif item.is_dir():
                    rmtree(item)  # Removes directories recursively
            logger.info("Directory '%s' has been deleted successfully.", dir_path)
        except FileNotFoundError:
            logger.info("Directory '%s' not found.", dir_path)
        except PermissionError as e:
            raise ValueError(
                f"Permission denied: Unable to delete '{dir_path}'."
            ) from e


def _align_obs_columns(experiment_df, anndata_df):
    def _create_na_columns(df_to_write, columns, dtype_df):
        """
        Create new columns in a dataframe with NA's as content while
        making sure that the datatype of the newly created columns match
        with those from another dataframe.
        """
        for col in columns:
            dtype = dtype_df[col].dtype
            na = pd.NA if pd.api.types.is_extension_array_dtype(dtype) else np.nan
            df_to_write.loc[:, col] = pd.Series(
                na, index=df_to_write.index, dtype=dtype
            )
        return df_to_write

    logger.info("Making sure obs columns are aligned.")
    contents_columns = set(experiment_df.columns)
    logger.info("Already ingested .obs columns: %s", ", ".join(contents_columns))
    logger.info("Retreiving obs columns to be added.")

    obs_to_add_columns = set(anndata_df.columns)

    logger.info("Checking dtype of common columns")
    common_columns = obs_to_add_columns.intersection(contents_columns)
    common_dtypes = {}
    for column_name in common_columns:
        # Combine the two columns, but only keep the calculated dtype.
        # We'll update the data types in the two original frames and let tiledb do the joining.
        new_series = experiment_df[column_name].combine_first(anndata_df[column_name])
        result_dtype = new_series.dtype
        # Nullable floats are not supported by anndata...
        if pd.api.types.is_float_dtype(
            result_dtype
        ) and pd.api.types.is_extension_array_dtype(result_dtype):
            result_dtype = result_dtype.numpy_dtype
        common_dtypes[column_name] = result_dtype

    logger.info(
        "Updating the dtypes to comply with the following schema: %s", common_dtypes
    )
    experiment_df, anndata_df = (
        experiment_df.astype(common_dtypes),
        anndata_df.astype(common_dtypes),
    )

    logger.info("Will add the following columns: %s", ", ".join(obs_to_add_columns))
    missing_columns_in_old = obs_to_add_columns - contents_columns
    logger.info(
        "Adding the following columns to the schema: %s.",
        ", ".join(missing_columns_in_old),
    )
    experiment_df = _create_na_columns(
        experiment_df, missing_columns_in_old, anndata_df
    )

    logger.info("Adjusting .obs succeeded!")
    missing_columns_in_new = (
        contents_columns
        - obs_to_add_columns
        - {"soma_joinid", par["obs_index_name_output"]}
    )
    logger.info(
        "Adding the following columns to modality to be ingested: %s",
        ", ".join(missing_columns_in_new),
    )
    anndata_df = _create_na_columns(anndata_df, missing_columns_in_new, experiment_df)
    return experiment_df, anndata_df


def main(par):
    logger.info(f"Component {meta['name']} started.")
    par["input"], par["tiledb_dir"] = Path(par["input"]), Path(par["tiledb_dir"])

    # Handle case where input is a mudat file
    with open(par["input"], "rb") as open_input_file:
        # It is possible to create a MuData compatible h5 file without using
        # MuData, and those could not have "MuData" as the first bytes.
        # But that is really an edge case and this check should hold.
        is_mudata = open_input_file.read(6) == b"MuData"

    _check_input_args(par, par["input"], is_mudata)

    # Reference layers other than "X" refer as "layer/<layer_name>"
    for arg_name in LAYER_ARGUMENTS:
        par[arg_name] = _handle_layer_location(par[arg_name])

    rna_input = par["input"]
    # Create a temporary file that holds the converted RNA h5ad
    if is_mudata:
        if not par["rna_modality"]:
            raise ValueError(
                "'rna_modality' argument must be set if the input is a MuData file."
            )
        # If the input is a MuData file, take the RNA modality and use that as input instead
        rna_input = _h5mu_to_h5ad(par["input"], par["rna_modality"])

    # Put a copy of the RNA input in the output file and convert it in-place
    rna_converted = _copy_anndata_and_convert_inplace(
        rna_input, _get_rna_conversion_specification(par)
    )

    if par["prot_modality"]:
        prot_input = _h5mu_to_h5ad(par["input"], par["prot_modality"])
        prot_converted = _copy_anndata_and_convert_inplace(
            prot_input, _get_prot_conversion_specification(par)
        )

        logger.info(
            "Making sure that obs columns are shared between the two modalities."
        )
        contents, obs_to_add = _read_obs(rna_converted), _read_obs(prot_converted)
        contents, obs_to_add = _align_obs_columns(contents, obs_to_add)
        _write_obs(prot_converted, obs_to_add)
        _write_obs(rna_converted, contents)

    tiledbsoma.logging.info()
    _remove_directory(par["tiledb_dir"])

    logger.info("Ingesting RNA modality")
    rna_from_h5ad_args = {
        "experiment_uri": str(par["tiledb_dir"]),
        "input_path": rna_converted,
        "measurement_name": par["rna_modality_output"],
        "uns_keys": [],
        "X_layer_name": "raw",
    }
    func_to_call = tiledbsoma.io.from_h5ad
    _log_arguments(func_to_call, rna_from_h5ad_args)
    func_to_call(**rna_from_h5ad_args)

    _write_varm_obsm_indices(par["tiledb_dir"], "rna", multidim_column_indices)

    logger.info("Done ingesting RNA modality")

    if par["prot_modality"]:
        logger.info("Ingesting protein modality")

        logger.info("Initalizing prot modality schema.")
        # Note the 'schema_only'
        prot_from_h5ad_schema_args = {
            "experiment_uri": str(par["tiledb_dir"]),
            "input_path": prot_converted,
            "measurement_name": par["prot_modality_output"],
            "uns_keys": [],
            "ingest_mode": "schema_only",
            "X_layer_name": "raw",
        }
        func_to_call = tiledbsoma.io.from_h5ad
        _log_arguments(func_to_call, prot_from_h5ad_schema_args)
        func_to_call(**prot_from_h5ad_schema_args)

        logger.info("Registering prot modality anndata")
        # Register the second anndata object in the protein measurement
        register_prot_args = {
            "experiment_uri": str(par["tiledb_dir"]),
            "h5ad_file_names": [prot_converted],
            "measurement_name": par["prot_modality_output"],
            "obs_field_name": par["obs_index_name_output"],
            "var_field_name": par["prot_var_index_name_output"],
            "append_obsm_varm": True,
        }
        func_to_call = tiledbsoma.io.register_h5ads
        _log_arguments(func_to_call, register_prot_args)
        rd = func_to_call(**register_prot_args)

        rd.prepare_experiment(str(par["tiledb_dir"]))

        # Ingest the second anndata object into the protein measurement
        prot_write_args = {
            "experiment_uri": str(par["tiledb_dir"]),
            "input_path": prot_converted,
            "measurement_name": par["prot_modality_output"],
            "registration_mapping": rd,
            "ingest_mode": "write",
            "uns_keys": [],
            "X_layer_name": "raw",
        }
        func_to_call = tiledbsoma.io.from_h5ad
        _log_arguments(func_to_call, prot_write_args)
        func_to_call(**prot_write_args)

        _write_varm_obsm_indices(par["tiledb_dir"], "prot", multidim_column_indices)

    logger.info("Finished!")


if __name__ == "__main__":
    main(par)
