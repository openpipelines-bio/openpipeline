from __future__ import annotations
import sys
import anndata
import mudata as mu
import pandas as pd
import numpy as np
from collections.abc import Iterable
from multiprocessing import Pool
from pathlib import Path
from h5py import File as H5File
from typing import Literal
import shutil

### VIASH START
par = {
    "input": [
        "resources_test/concat_test_data/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix_subset_unique_obs.h5mu",
        "resources_test/concat_test_data/human_brain_3k_filtered_feature_bc_matrix_subset_unique_obs.h5mu",
    ],
    "output": "foo.h5mu",
    "input_id": ["mouse", "human"],
    "obsp_keys": [],
    "other_axis_mode": "move",
    "output_compression": "gzip",
    "uns_merge_mode": "make_unique",
}
meta = {"cpus": 10, "resources_dir": "resources_test/"}
### VIASH END

sys.path.append(meta["resources_dir"])
from compress_h5mu import compress_h5mu
from setup_logger import setup_logger

logger = setup_logger()

# 'move' is handled after concatenation, start from the 'unique' merge strategy
OTHER_AXIS_CONCAT_MODES = {"move": "unique"}
# 'make_unique' is handled after concatenation, start from not merging any keys
UNS_CONCAT_MODES = {"make_unique": None}


def nunique(row):
    unique = pd.unique(row)
    unique_without_na = pd.core.dtypes.missing.remove_na_arraylike(unique)
    return len(unique_without_na) > 1


def any_row_contains_duplicate_values(n_processes: int, frame: pd.DataFrame) -> bool:
    """
    Check if any row contains duplicate values, that are not NA.
    """
    numpy_array = frame.to_numpy()
    with Pool(n_processes) as pool:
        is_duplicated = pool.map(nunique, iter(numpy_array))
    return any(is_duplicated)


def concatenate_matrices(
    n_processes: int, matrices: dict[str, pd.DataFrame], align_to: pd.Index
) -> tuple[
    dict[str, pd.DataFrame], pd.DataFrame | None, dict[str, pd.core.dtypes.dtypes.Dtype]
]:
    """
    Merge matrices by combining columns that have the same name.
    Columns that contain conflicting values (e.i. the columns have different values),
    are not merged, but instead moved to a new dataframe.
    """
    column_names = set(column_name for var in matrices.values() for column_name in var)
    logger.debug("Trying to concatenate columns: %s.", ",".join(column_names))
    if not column_names:
        return {}, pd.DataFrame(index=align_to)
    conflicts, concatenated_matrix = split_conflicts_and_concatenated_columns(
        n_processes, matrices, column_names, align_to
    )
    concatenated_matrix = cast_to_writeable_dtype(concatenated_matrix)
    conflicts = {
        conflict_name: cast_to_writeable_dtype(conflict_df)
        for conflict_name, conflict_df in conflicts.items()
    }
    return conflicts, concatenated_matrix


def get_first_non_na_value_vector(df):
    numpy_arr = df.to_numpy()
    n_rows, n_cols = numpy_arr.shape
    col_index = pd.isna(numpy_arr).argmin(axis=1)
    flat_index = n_cols * np.arange(n_rows) + col_index
    return pd.Series(numpy_arr.ravel()[flat_index], index=df.index, name=df.columns[0])


def make_uns_keys_unique(mod_data, concatenated_data):
    """
    Check if the uns keys across samples are unique before adding them
    to the final concatenated object. If a conflict occurs between the samples,
    add the sample ID to make the key unique again.
    """
    all_uns_keys = {}
    for sample_id, mod in mod_data.items():
        for uns_key, _ in mod.uns.items():
            all_uns_keys.setdefault(uns_key, []).append(sample_id)
    for uns_key, samples_ids in all_uns_keys.items():
        assert samples_ids
        if len(samples_ids) == 1:
            sample_id = samples_ids[0]
            concatenated_data.uns[uns_key] = mod_data[sample_id].uns[uns_key]
        else:
            for sample_id in samples_ids:
                concatenated_data.uns[f"{sample_id}_{uns_key}"] = mod_data[
                    sample_id
                ].uns[uns_key]
    return concatenated_data


def split_conflicts_and_concatenated_columns(
    n_processes: int,
    matrices: dict[str, pd.DataFrame],
    column_names: Iterable[str],
    align_to: pd.Index,
) -> tuple[dict[str, pd.DataFrame], pd.DataFrame]:
    """
    Retrieve columns with the same name from a list of dataframes which are
    identical across all the frames (ignoring NA values).
    Columns which are not the same are regarded as 'conflicts',
    which are stored in seperate dataframes, one per columns
    with the same name that store conflicting values.
    """
    conflicts = {}
    concatenated_matrix = []
    for column_name in column_names:
        columns = {
            input_id: var[column_name]
            for input_id, var in matrices.items()
            if column_name in var
        }
        assert columns, "Some columns should have been found."
        concatenated_columns = pd.concat(
            columns.values(), axis=1, join="outer", sort=False
        )
        if any_row_contains_duplicate_values(n_processes, concatenated_columns):
            concatenated_columns.columns = (
                columns.keys()
            )  # Use the sample id as column name
            concatenated_columns = concatenated_columns.reindex(align_to, copy=False)
            conflicts[f"conflict_{column_name}"] = concatenated_columns
        else:
            unique_values = get_first_non_na_value_vector(concatenated_columns)
            concatenated_matrix.append(unique_values)
    if not concatenated_matrix:
        return conflicts, pd.DataFrame(index=align_to)
    concatenated_matrix = pd.concat(
        concatenated_matrix, join="outer", axis=1, sort=False
    )
    concatenated_matrix = concatenated_matrix.reindex(align_to, copy=False)
    return conflicts, concatenated_matrix


def cast_to_writeable_dtype(result: pd.DataFrame) -> pd.DataFrame:
    """
    Cast the dataframe to dtypes that can be written by mudata.
    """
    # dtype inferral workfs better with np.nan
    result = result.replace({pd.NA: np.nan})

    # MuData supports nullable booleans and ints
    # ie. `IntegerArray` and `BooleanArray`
    result = result.convert_dtypes(
        infer_objects=True,
        convert_integer=True,
        convert_string=False,
        convert_boolean=True,
        convert_floating=False,
    )

    # Convert leftover 'object' columns (strings or mixed types) to the 'str' dtype,
    # which keeps NA values and can be written by anndata
    object_cols = result.columns[result.dtypes == "object"]
    result[object_cols] = result[object_cols].astype("str")
    return result


def split_conflicts_modalities(
    n_processes: int, samples: dict[str, anndata.AnnData], output: anndata.AnnData
) -> anndata.AnnData:
    """
    Merge .var and .obs matrices of the anndata objects. Columns are merged
    when the values (excl NA) are the same in each of the matrices.
    Conflicting columns are moved to a separate dataframe (one dataframe for each column,
    containing all the corresponding column from each sample).
    """
    matrices_to_parse = ("var", "obs")
    for matrix_name in matrices_to_parse:
        matrices = {
            sample_id: getattr(sample, matrix_name)
            for sample_id, sample in samples.items()
        }
        output_index = getattr(output, matrix_name).index
        conflicts, concatenated_matrix = concatenate_matrices(
            n_processes, matrices, output_index
        )
        if concatenated_matrix.empty:
            concatenated_matrix.index = output_index

        # Even though we did not touch the varm and obsm matrices that were already present,
        # the joining of observations might have caused a dtype change in these matrices as well
        # so these also need to be casted to a writable dtype...
        for multidim_name, multidim_data in getattr(output, f"{matrix_name}m").items():
            new_data = (
                cast_to_writeable_dtype(multidim_data)
                if isinstance(multidim_data, pd.DataFrame)
                else multidim_data
            )
            getattr(output, f"{matrix_name}m")[multidim_name] = new_data

        # Write the conflicts to the output
        for conflict_name, conflict_data in conflicts.items():
            getattr(output, f"{matrix_name}m")[conflict_name] = conflict_data

        # Set other annotation matrices in the output
        setattr(output, matrix_name, concatenated_matrix)

    return output


def concatenate_modality(
    n_processes: int,
    mod: str,
    input_files: Iterable[str | Path],
    other_axis_mode: str,
    uns_merge_mode: str,
    input_ids: tuple[str],
) -> anndata.AnnData:
    other_axis_mode_to_apply = OTHER_AXIS_CONCAT_MODES.get(
        other_axis_mode, other_axis_mode
    )
    uns_merge_mode_to_apply = UNS_CONCAT_MODES.get(uns_merge_mode, uns_merge_mode)

    mod_data = {}
    for input_id, input_file in zip(input_ids, input_files):
        try:
            data = mu.read_h5ad(input_file, mod=mod)

            # Remove obsp keys that are not in par["obsp_keys"]
            obsp_keys_to_keep = par.get("obsp_keys") or []
            obsp_keys_to_remove = set(data.obsp.keys()) - set(obsp_keys_to_keep)
            for key in obsp_keys_to_remove:
                try:
                    del data.obsp[key]
                except KeyError:
                    pass

            mod_data[input_id] = data

        except KeyError as e:  # Modality does not exist for this sample, skip it
            if (
                f"Unable to synchronously open object (object '{mod}' doesn't exist)"
                not in str(e)
            ):
                raise e
            pass

    if not mod_data:
        return anndata.AnnData()

    concatenated_data = anndata.concat(
        mod_data.values(),
        join="outer",
        pairwise=True if par["obsp_keys"] else False,
        merge=other_axis_mode_to_apply,
        uns_merge=uns_merge_mode_to_apply,
    )

    if other_axis_mode == "move":
        concatenated_data = split_conflicts_modalities(
            n_processes, mod_data, concatenated_data
        )

    if uns_merge_mode == "make_unique":
        concatenated_data = make_uns_keys_unique(mod_data, concatenated_data)

    return concatenated_data


def get_modality_obs_names(input_file: Path, modalities: Iterable[str]) -> pd.Index:
    """
    Get the observation names of the selected modalities of a sample.
    """
    obs_names = pd.Index([])
    with H5File(input_file, "r") as input_h5:
        for mod in modalities:
            if mod not in input_h5["mod"]:
                continue
            obs_group = input_h5["mod"][mod]["obs"]
            index = anndata.io.read_elem(obs_group[obs_group.attrs["_index"]])
            obs_names = obs_names.union(pd.Index(index))
    return obs_names


def split_global_columns(
    column_names: Iterable[str],
    modality_columns: dict[str, set[str]],
    excluded_modalities: set[str],
) -> tuple[set[str], list[str]]:
    """
    Split the column names of a global annotation frame into the names of the
    modality columns they were pulled from, and columns that only exist globally.

    A column matches a modality column by its base name: either 'mod:column'
    (split on the first colon, as mudata adds a single modality prefix) or 'column'.
    Columns that refer to a modality that is not being output are dropped.
    """
    all_modality_columns = set().union(*modality_columns.values())
    to_pull, global_only = set(), []
    for column_name in column_names:
        prefix, sep, base_name = column_name.partition(":")
        if sep and prefix in excluded_modalities:
            logger.info(
                "Dropping global column '%s': modality '%s' is not being output.",
                column_name,
                prefix,
            )
        elif sep and base_name in modality_columns.get(prefix, ()):
            to_pull.add(base_name)
        elif column_name in all_modality_columns:
            to_pull.add(column_name)
        else:
            global_only.append(column_name)
    return to_pull, global_only


def concatenate_global_annotation(
    n_processes: int,
    mdata: mu.MuData,
    axis: Literal["obs", "var"],
    input_files: list[Path],
    input_ids: tuple[str],
    sample_obs_names: dict[str, pd.Index],
    other_axis_mode: str,
) -> None:
    """
    Add the global annotation columns (.obs or .var) from the samples to mdata.
    Columns that were pulled from a modality are pulled again from the concatenated
    modalities, using the current modality values and the mudata naming
    ('mod:column', or 'column' for .var columns that are present in all modalities).
    Columns that only exist globally are concatenated across the samples.
    """
    modality_columns = {
        mod: set(getattr(mod_data, axis).columns) for mod, mod_data in mdata.mod.items()
    }
    columns_to_pull = set()
    samples = {}
    for input_id, input_file in zip(input_ids, input_files):
        with H5File(input_file, "r") as input_h5:
            excluded_modalities = set(input_h5["mod"].keys()) - set(mdata.mod.keys())
            frame = (
                anndata.io.read_elem(input_h5[axis])
                if axis in input_h5
                else pd.DataFrame()
            )
        to_pull, global_only = split_global_columns(
            frame.columns, modality_columns, excluded_modalities
        )
        columns_to_pull |= to_pull
        frame = frame[global_only]
        if axis == "obs":
            # Only keep the observations from the modalities that are being output
            frame = frame[frame.index.isin(sample_obs_names[input_id])]
        samples[input_id] = anndata.AnnData(**{axis: frame})

    if columns_to_pull:
        getattr(mdata, f"pull_{axis}")(columns=sorted(columns_to_pull))

    concatenated_data = anndata.concat(
        samples.values(),
        join="outer",
        merge=OTHER_AXIS_CONCAT_MODES.get(other_axis_mode, other_axis_mode),
    )
    if other_axis_mode == "move":
        concatenated_data = split_conflicts_modalities(
            n_processes, samples, concatenated_data
        )

    names = getattr(mdata, f"{axis}_names")
    global_only = cast_to_writeable_dtype(
        getattr(concatenated_data, axis).reindex(names)
    )
    setattr(mdata, axis, pd.concat([getattr(mdata, axis), global_only], axis=1))
    for key, value in getattr(concatenated_data, f"{axis}m").items():
        getattr(mdata, f"{axis}m")[key] = cast_to_writeable_dtype(value.reindex(names))


def concatenate_global_uns(
    input_files: list[Path], input_ids: tuple[str], uns_merge_mode: str
) -> dict:
    """
    Merge the global .uns of the samples.
    """
    samples = {}
    for input_id, input_file in zip(input_ids, input_files):
        with H5File(input_file, "r") as input_h5:
            uns = anndata.io.read_elem(input_h5["uns"]) if "uns" in input_h5 else {}
        if uns:
            samples[input_id] = anndata.AnnData(uns=uns)
    if not samples:
        return {}

    concatenated_data = anndata.concat(
        samples.values(),
        uns_merge=UNS_CONCAT_MODES.get(uns_merge_mode, uns_merge_mode),
    )
    if uns_merge_mode == "make_unique":
        concatenated_data = make_uns_keys_unique(samples, concatenated_data)
    return dict(concatenated_data.uns)


def concatenate_global_data(
    n_processes: int,
    modalities: list[str],
    input_files: list[Path],
    input_ids: tuple[str],
    sample_obs_names: dict[str, pd.Index],
    other_axis_mode: str,
    uns_merge_mode: str,
    output_file: Path,
) -> None:
    """
    Add the global .obs, .var and .uns from the samples to the output.
    """
    logger.info("Concatenating global data.")
    # Only the annotation frames are needed to determine the global obs and var names
    with H5File(output_file, "r") as output_h5:
        mdata = mu.MuData(
            {
                mod: anndata.AnnData(
                    obs=anndata.io.read_elem(output_h5["mod"][mod]["obs"]),
                    var=anndata.io.read_elem(output_h5["mod"][mod]["var"]),
                )
                for mod in modalities
            }
        )

    with H5File(output_file, "r+") as output_h5:
        for axis in ("obs", "var"):
            concatenate_global_annotation(
                n_processes,
                mdata,
                axis,
                input_files,
                input_ids,
                sample_obs_names,
                other_axis_mode,
            )
            # Write the global elements the same way mudata does
            anndata.io.write_elem(
                output_h5,
                axis,
                mdata.strings_to_categoricals(getattr(mdata, axis).copy()),
            )
            for element in (f"{axis}m", f"{axis}map"):
                anndata.io.write_elem(output_h5, element, dict(getattr(mdata, element)))
        anndata.io.write_elem(
            output_h5,
            "uns",
            concatenate_global_uns(input_files, input_ids, uns_merge_mode),
        )


def concatenate_modalities(
    n_processes: int,
    modalities: list[str],
    input_files: Path | str,
    other_axis_mode: str,
    uns_merge_mode: str,
    output_file: Path | str,
    compression: Literal["gzip"] | Literal["lzf"],
    input_ids: tuple[str] | None = None,
) -> None:
    """
    Join the modalities together into a single multimodal sample.
    """
    logger.info("Concatenating samples.")
    output_file, input_files = (
        Path(output_file),
        [Path(input_file) for input_file in input_files],
    )
    # The same observation can be in different modalities of different samples,
    # which is not detected when checking each modality separately.
    sample_obs_names = {
        input_id: get_modality_obs_names(input_file, modalities)
        for input_id, input_file in zip(input_ids, input_files)
    }
    if not pd.Index([]).append(list(sample_obs_names.values())).is_unique:
        raise ValueError("Observations are not unique across samples.")

    output_file_uncompressed = output_file.with_name(
        output_file.stem + "_uncompressed.h5mu"
    )
    output_file_uncompressed.touch()
    # Create empty mudata file
    mdata = mu.MuData({modality: anndata.AnnData() for modality in modalities})
    mdata.write(output_file_uncompressed, compression=compression)

    for mod_name in modalities:
        new_mod = concatenate_modality(
            n_processes,
            mod_name,
            input_files,
            other_axis_mode,
            uns_merge_mode,
            input_ids,
        )
        logger.info(
            "Writing out modality '%s' to '%s' with compression '%s'.",
            mod_name,
            output_file_uncompressed,
            compression,
        )
        mu.write_h5ad(output_file_uncompressed, data=new_mod, mod=mod_name)

    concatenate_global_data(
        n_processes,
        modalities,
        input_files,
        input_ids,
        sample_obs_names,
        other_axis_mode,
        uns_merge_mode,
        output_file_uncompressed,
    )

    if compression:
        compress_h5mu(output_file_uncompressed, output_file, compression=compression)
        output_file_uncompressed.unlink()
    else:
        shutil.move(output_file_uncompressed, output_file)

    logger.info("Concatenation successful.")


def main() -> None:
    # Get a list of all possible modalities
    mods = set()
    for path in par["input"]:
        try:
            with H5File(path, "r") as f_root:
                mods = mods | set(f_root["mod"].keys())
        except OSError:
            raise OSError(f"Failed to load {path}. Is it a valid h5 file?")

    input_ids = None
    if par["input_id"]:
        input_ids: tuple[str] = tuple(i.strip() for i in par["input_id"])
        if len(input_ids) != len(par["input"]):
            raise ValueError(
                "The number of sample names must match the number of sample files."
            )

        if len(set(input_ids)) != len(input_ids):
            raise ValueError("The sample names should be unique.")

    logger.info("\nConcatenating data from paths:\n\t%s", "\n\t".join(par["input"]))

    if par["other_axis_mode"] == "move" and not input_ids:
        raise ValueError("--mode 'move' requires --input_ids.")

    n_processes = meta["cpus"] if meta["cpus"] else 1

    if par["modality"]:
        par["modality"] = set(par["modality"])
        if not par["modality"].issubset(mods):
            mods_joined, input_mods_joined = ", ".join(mods), ", ".join(par["modality"])
            raise ValueError(
                f"One of the modalities provided ({input_mods_joined}) is not available in the input data {mods_joined}"
            )
        mods = par["modality"]

    concatenate_modalities(
        n_processes,
        list(mods),
        par["input"],
        par["other_axis_mode"],
        par["uns_merge_mode"],
        par["output"],
        par["output_compression"],
        input_ids=input_ids,
    )


if __name__ == "__main__":
    main()
