from __future__ import annotations
import sys
import anndata
import mudata as mu
import pandas as pd
import numpy as np
from collections.abc import Iterable
from multiprocessing import Pool
from pathlib import Path
import h5py

### VIASH START
par = {
    "input": ["resources_test/concat_test_data/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix_subset_unique_obs.h5mu",
              "resources_test/concat_test_data/human_brain_3k_filtered_feature_bc_matrix_subset_unique_obs.h5mu"],
    "output": "foo.h5mu",
    "input_id": ["mouse", "human"],
    "compression": "gzip",
    "other_axis_mode": "move",
    "output_compression": "gzip"
}
meta = {
    "cpus": 10
}
### VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
logger = setup_logger()

def indexes_unique(indices: Iterable[pd.Index]) -> bool:
    combined_indices = indices[0].append(indices[1:])
    return combined_indices.is_unique

def check_observations_unique(samples: Iterable[anndata.AnnData]) -> None:
    observation_ids = [sample.obs.index for sample in samples]
    if not indexes_unique(observation_ids):
        raise ValueError("Observations are not unique across samples.")


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

def concatenate_matrices(n_processes: int, input_ids: tuple[str], matrices: Iterable[pd.DataFrame]) \
    -> tuple[dict[str, pd.DataFrame], pd.DataFrame | None, dict[str, pd.core.dtypes.dtypes.Dtype]]:
    """
    Merge matrices by combining columns that have the same name.
    Columns that contain conflicting values (e.i. the columns have different values),
    are not merged, but instead moved to a new dataframe.
    """
    column_names = set(column_name for var in matrices for column_name in var)
    logger.debug('Trying to concatenate columns: %s.', ",".join(column_names))
    if not column_names:
        return {}, None
    conflicts, concatenated_matrix = \
        split_conflicts_and_concatenated_columns(n_processes,
                                                 input_ids,
                                                 matrices,
                                                 column_names)
    concatenated_matrix = cast_to_writeable_dtype(concatenated_matrix)
    conflicts = {conflict_name: cast_to_writeable_dtype(conflict_df) 
                 for conflict_name, conflict_df in conflicts.items()}
    return conflicts, concatenated_matrix

def get_first_non_na_value_vector(df):
    numpy_arr = df.to_numpy()
    n_rows, n_cols = numpy_arr.shape
    col_index = pd.isna(numpy_arr).argmin(axis=1)
    flat_index = n_cols * np.arange(n_rows) + col_index
    return pd.Series(numpy_arr.ravel()[flat_index], index=df.index, name=df.columns[0])

def split_conflicts_and_concatenated_columns(n_processes: int,
                                             input_ids: tuple[str],
                                             matrices: Iterable[pd.DataFrame],
                                             column_names: Iterable[str]) -> \
                                            tuple[dict[str, pd.DataFrame], pd.DataFrame]:
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
        columns = [var[column_name] for var in matrices if column_name in var]
        assert columns, "Some columns should have been found."
        concatenated_columns = pd.concat(columns, axis=1, join="outer")
        if any_row_contains_duplicate_values(n_processes, concatenated_columns):
            concatenated_columns.columns = input_ids
            conflicts[f'conflict_{column_name}'] = concatenated_columns
        else:
            unique_values = get_first_non_na_value_vector(concatenated_columns)
            # concatenated_columns.fillna(method='bfill', axis=1).iloc[:, 0]
            concatenated_matrix.append(unique_values)
    if concatenated_matrix:
        concatenated_matrix = pd.concat(concatenated_matrix, join="outer", axis=1)
    else:
        concatenated_matrix = pd.DataFrame()

    return conflicts, concatenated_matrix

def cast_to_writeable_dtype(result: pd.DataFrame) -> pd.DataFrame:
    """
    Cast the dataframe to dtypes that can be written by mudata.
    """    
    # dtype inferral workfs better with np.nan
    result = result.replace({pd.NA: np.nan})

    # MuData supports nullable booleans and ints
    # ie. `IntegerArray` and `BooleanArray`
    result = result.convert_dtypes(infer_objects=True,
                                   convert_integer=True,
                                   convert_string=False,
                                   convert_boolean=True,
                                   convert_floating=False)
    
    # Convert leftover 'object' columns to string
    object_cols = result.select_dtypes(include='object').columns.values
    for obj_col in object_cols:
        result[obj_col].astype(str).astype('category')
    return result

def split_conflicts_modalities(n_processes: int, input_ids: tuple[str], modalities: Iterable[anndata.AnnData]) \
        -> tuple[dict[str, dict[str, pd.DataFrame]],  dict[str, pd.DataFrame | None]]:
        """
        Merge .var and .obs matrices of the anndata objects. Columns are merged
        when the values (excl NA) are the same in each of the matrices.
        Conflicting columns are moved to a separate dataframe (one dataframe for each column,
        containing all the corresponding column from each sample).
        """
        matrices_to_parse = ("var", "obs")
        concatenated_result = {}
        conflicts_result = {}
        for matrix_name in matrices_to_parse:
            matrices = [getattr(modality, matrix_name) for modality in modalities]
            conflicts, concatenated_matrix = concatenate_matrices(n_processes, input_ids, matrices)
            conflicts_result[f"{matrix_name}m"] = conflicts
            concatenated_result[matrix_name] = concatenated_matrix
        return conflicts_result, concatenated_result

def set_matrices(mod: mu.MuData,
                 new_matrices: dict[str, pd.DataFrame | None]) -> mu.MuData:
    """
    Add the calculated matrices to the mudata object. Ensure the correct datatypes
    for the matrices that are composed from the combination of the matrices
    from the different modalities.
    """
    for matrix_name, data in new_matrices.items():
        new_index = getattr(mod, matrix_name).index
        if data is None:
            data = pd.DataFrame(index=new_index)
        if data.index.empty:
            data.index = new_index
        setattr(mod, matrix_name, data)
    return mod


def set_conflicts(mod: anndata.AnnData,
                  conflicts: dict[str, dict[str, pd.DataFrame]]) -> mu.MuData:
    """
    Store dataframes containing the conflicting columns in .obsm,
    one key per column name from the original data.
    """
    mutlidim_to_singledim = {
        'varm': 'var',
        'obsm': 'obs'
    }
    for conflict_matrix_name, conflict in conflicts.items():
        for conflict_name, conflict_data in conflict.items():
            singledim_name = mutlidim_to_singledim[conflict_matrix_name]
            singledim_index = getattr(mod, singledim_name).index
            getattr(mod, conflict_matrix_name)[conflict_name] = conflict_data.reindex(singledim_index)
    return mod

def concatenate_modality(n_processes: int, mod: str, input_files: Iterable[str | Path], 
                         other_axis_mode: str, input_ids: tuple[str]) -> anndata.AnnData:
    
    concat_modes = {
        "move": None,
    }
    other_axis_mode_to_apply = concat_modes.get(other_axis_mode, other_axis_mode)

    mod_data = []
    for input_file in input_files:
        try:
            mod_data.append(mu.read_h5ad(input_file, mod=mod))
        except KeyError as e: # Modality does not exist for this sample, skip it
            if f"Unable to open object '{mod}' doesn't exist" not in str(e):
                raise e
            pass
    check_observations_unique(mod_data)

    concatenated_data = anndata.concat(mod_data, join='outer', merge=other_axis_mode_to_apply)

    if other_axis_mode == "move":
        conflicts, new_matrices = split_conflicts_modalities(n_processes, input_ids, mod_data)
        concatenated_data = set_conflicts(concatenated_data, conflicts)
        concatenated_data = set_matrices(concatenated_data, new_matrices)
    return concatenated_data

def concatenate_modalities(n_processes: int, modalities: set[str], input_files: Path | str,
                           other_axis_mode: str, input_ids: tuple[str] | None = None) -> mu.MuData:
    """
    Join the modalities together into a single multimodal sample.
    """
    logger.info('Concatenating samples.')
    if other_axis_mode == "move" and not input_ids:
        raise ValueError("--mode 'move' requires --input_ids.")
    new_mods = {}
    for mod_name in modalities:
        new_mods[mod_name] = concatenate_modality(n_processes, mod_name, input_files, other_axis_mode, input_ids)
    concatenated_data = mu.MuData(new_mods)

    # After setting the matrices (.e.g. .mod['rna'].var) for each of the modalities
    # the 'global' matrices must also be updated. This is done by mudata automatically
    # by calling mudata.update() before writing. However, we need to make sure that the
    # dtypes of these global matrices are also correct for writing..
    for global_matrix_name in ("var", "obs"):
        matrix = getattr(concatenated_data, global_matrix_name)
        setattr(concatenated_data, global_matrix_name, cast_to_writeable_dtype(matrix))
    logger.info("Concatenation successful.")
    return concatenated_data


def main() -> None:
    # Get a list of all possible modalities
    mods = set()
    for path in par["input"]:
        try:
            with h5py.File(path, 'r') as f_root:
                mods = mods | set(f_root["mod"].keys())
        except OSError:
            raise OSError(f"Failed to load {path}. Is it a valid h5 file?")

    input_ids = None
    if par["input_id"]:
        input_ids: tuple[str] = tuple(i.strip() for i in par["input_id"])
        if len(input_ids) != len(par["input"]):
            raise ValueError("The number of sample names must match the number of sample files.")

        if len(set(input_ids)) != len(input_ids):
            raise ValueError("The sample names should be unique.")

    logger.info("\nConcatenating data from paths:\n\t%s",
                "\n\t".join(par["input"]))

    n_processes = meta["cpus"] if meta["cpus"] else 1
    concatenated_samples = concatenate_modalities(n_processes,
                                                  mods,
                                                  par["input"],
                                                  par["other_axis_mode"],
                                                  input_ids=input_ids)
    logger.info("Writing out data to '%s' with compression '%s'.",
                par["output"], par["output_compression"])
    concatenated_samples.write_h5mu(par["output"], compression=par["output_compression"])


if __name__ == "__main__":
    main()
