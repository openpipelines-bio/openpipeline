from __future__ import annotations
import logging
import anndata
import muon as mu
from sys import stdout
import pandas as pd
import numpy as np
from collections.abc import Iterable
from multiprocessing import Pool

### VIASH START
par = {
    "input": ["resources_test/concat/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix_subset.h5mu",
              "resources_test/concat/human_brain_3k_filtered_feature_bc_matrix_subset.h5mu"],
    "output": "foo.h5mu",
    "sample_names": ["mouse", "human"],
    "obs_sample_name": "sample_id",
    "compression": "gzip",
    "other_axis_mode": "move"
}
meta = {
    "cpus": 10

}
### VIASH END

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

def add_sample_names(sample_ids: tuple[str], samples: Iterable[mu.MuData], obs_key_sample_name: str) -> None:
    """
    Add sample names to the observations for each sample.
    """
    for (sample_id, sample) in zip(sample_ids, samples):
        if obs_key_sample_name in sample.obs_keys():
            logger.info(f'Column .obs["{obs_key_sample_name}"] already exists in sample "{sample_id}". Overriding the value for this column.')
            sample.obs = sample.obs.drop(obs_key_sample_name, axis=1)
        for modality in sample.mod.values():
            modality.obs[obs_key_sample_name] = sample_id
        sample.update()


def make_observation_keys_unique(sample_ids: tuple[str], samples: Iterable[mu.MuData]) -> None:
    """
    Make the observation keys unique across all samples. At input,
    the observation keys are unique within a sample. By adding the sample name
    (unique for a sample) to each observation key, the observation key is made
    unique across all samples as well.
    """
    logger.info('Making observation keys unique across all samples.')
    for sample_id, sample in zip(sample_ids, samples):
        sample.obs.index = f"{sample_id}_" + sample.obs.index
        make_observation_keys_unique_per_mod(sample_id, sample)


def make_observation_keys_unique_per_mod(sample_id: str, sample: mu.MuData) -> None:
    """
    Updating MuData.obs_names is not allowed (it is read-only).
    So the observation keys for each modality has to be updated manually.
    """
    for mod in sample.mod.values():
        mod.obs_names = f"{sample_id}_" + mod.obs_names


def group_modalities(samples: Iterable[anndata.AnnData]) -> dict[str, anndata.AnnData]:
    """
    Split up the modalities of all samples and group them per modality.
    """
    mods = {}
    for sample in samples:
        for mod_name, mod in sample.mod.items():
            mods.setdefault(mod_name, []).append(mod)

    if len(set(len(mod) for mod in mods.values())) != 1:
        logger.warning("One or more samples seem to have a different number of modalities.")

    logger.info("Successfully sorted modalities for the different samples.")
    return mods

def nunique(row):
    unique = pd.unique(row)
    unique_without_na = pd.core.dtypes.missing.remove_na_arraylike(unique)
    return len(unique_without_na)

def any_row_contains_duplicate_values(n_processes: int, frame: pd.DataFrame) -> bool:
    """
    Check if any row contains duplicate values, that are not NA.
    """
    numpy_array = frame.to_numpy()
    with Pool(n_processes) as pool:
        number_of_unique = pool.map(nunique, iter(numpy_array))
    number_of_unique = pd.Series(number_of_unique, index=frame.index)
    non_na_counts = frame.count(axis=1)
    is_duplicated = (number_of_unique - non_na_counts) != 0
    return is_duplicated.any()

def concatenate_matrices(n_processes: int, sample_ids: tuple[str], matrices: Iterable[pd.DataFrame]) \
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
                                                 sample_ids,
                                                 matrices,
                                                 column_names)
    original_dtypes = get_original_dtypes(matrices, column_names)
    concatenated_matrix = set_dtypes_concatenated_columns(original_dtypes, concatenated_matrix)
    conflicts = set_dtypes_conflicts(original_dtypes, conflicts)
    return conflicts, concatenated_matrix

def set_dtypes_conflicts(original_dtypes: dict[str, pd.core.dtypes.dtypes.Dtype],
                         conflicts: tuple[dict[str, pd.DataFrame], pd.DataFrame]) -> \
                        tuple[dict[str, pd.DataFrame], pd.DataFrame]:
    """
    Ensure the correct datatypes for the conflict dataframes.
    """
    conflicts_correct_dtypes = {}
    for conflict_name, confict_data in conflicts.items():
        original_dtype = original_dtypes[conflict_name.removeprefix('conflict_')]
        new_conflict_dtypes = {column: original_dtype for column in confict_data.columns}
        confict_data = cast_to_original_dtype(confict_data, new_conflict_dtypes)
        conflicts_correct_dtypes[conflict_name] = confict_data
    return conflicts_correct_dtypes

def set_dtypes_concatenated_columns(original_dtypes: dict[str, pd.core.dtypes.dtypes.Dtype],
                                    concatenated_matrix: pd.DataFrame) -> pd.DataFrame:
    """
    Ensure the correct datatypes for the concatenated columns that did not contain conflicts.
    """
    curr_concat_matrix_cols_dtypes = {col: dtype for col, dtype in original_dtypes.items() 
                                      if col in concatenated_matrix.columns}
    return cast_to_original_dtype(concatenated_matrix, curr_concat_matrix_cols_dtypes)

def get_first_non_na_value_vector(df):
    numpy_arr = df.to_numpy()
    n_rows, n_cols = numpy_arr.shape
    col_index = pd.isna(numpy_arr).argmin(axis=1)
    flat_index = n_cols * np.arange(n_rows) + col_index
    return pd.Series(numpy_arr.ravel()[flat_index], index=df.index, name=df.columns[0])

def split_conflicts_and_concatenated_columns(n_processes: int,
                                             sample_ids: tuple[str],
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
            concatenated_columns.columns = sample_ids
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

def cast_to_original_dtype(result: pd.DataFrame,
                           orignal_dtypes: dict[str, pd.core.dtypes.dtypes.Dtype]) -> pd.DataFrame:
    """
    Cast the dataframe to dtypes that can be written by mudata.
    """
    logger.debug('Trying to cast to "category" or keep original datatype.')
    for col_name, orig_dtype in orignal_dtypes.items():
        try:
            result = result.astype({col_name: "category"}, copy=True)
            result[col_name].cat.categories = result[col_name].cat.categories.astype(str)
        except (ValueError, TypeError):
            try:
                result = result.astype({col_name: orig_dtype}, copy=True)
            except (ValueError, TypeError):
                logger.warning("Could not keep datatype for column %s", col_name)
    return result


def get_original_dtypes(matrices: Iterable[pd.DataFrame], 
                        column_names: Iterable[str]) -> \
                        dict[str, pd.core.dtypes.dtypes.Dtype]:
    """
    Get the datatypes of columns in a list of dataframes.
    If a column occurs in more than 1 dataframe, includes the dtype of the column
    in the dataframe that comes first in the list.
    """
    dtypes = {}
    for col_name in column_names:
        for matrix in matrices:
            col = matrix.get(col_name, None)
            if col is not None and col_name not in dtypes:
                dtypes[col_name] = col.dtype
    return dtypes


def split_conflicts_modalities(n_processes: int, sample_ids: tuple[str], modalities: Iterable[anndata.AnnData]) \
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
            conflicts, concatenated_matrix = concatenate_matrices(n_processes, sample_ids, matrices)
            conflicts_result[f"{matrix_name}m"] = conflicts
            concatenated_result[matrix_name] = concatenated_matrix
        return conflicts_result, concatenated_result

def set_matrices(concatenated_data: mu.MuData, 
                 mod_name: str,
                 new_matrices: dict[str, pd.DataFrame | None]) -> mu.MuData:
    """
    Add the calculated matrices to the mudata object. Ensure the correct datatypes
    for the matrices that are composed from the combination of the matrices
    from the different modalities.
    """
    mod = concatenated_data.mod[mod_name]
    original_dtypes_global_matrices = {
        global_matrix_name: getattr(concatenated_data, global_matrix_name).dtypes
        for global_matrix_name 
        in new_matrices.keys()
    }
    for matrix_name, data in new_matrices.items():
        new_index = getattr(mod, matrix_name).index
        if data is None:
            data = pd.DataFrame(index=new_index)
        if data.index.empty:
            data.index = new_index
        setattr(mod, matrix_name, data)
    # After setting the matrices (.e.g. .mod['rna'].var) for each of the modalities
    # the 'global' matrices must also be updated. This is done by mudata automatically
    # by calling mudata.update() before writing. However, we need to make sure that the
    # dtypes of these global matrices are also correct for writing..
    for global_matrix_name, dtypes in original_dtypes_global_matrices.items():
        matrix = getattr(concatenated_data, global_matrix_name)
        setattr(concatenated_data, global_matrix_name, cast_to_original_dtype(matrix, dtypes))
    return concatenated_data

def set_conflicts(concatenated_data: mu.MuData,
                  mod_name: str,
                  conflicts: dict[str, dict[str, pd.DataFrame]]) -> mu.MuData:
    """
    Store dataframes containing the conflicting columns in .obsm,
    one key per column name from the original data.
    """
    mod = concatenated_data.mod[mod_name]
    mutlidim_to_singledim = {
        'varm': 'var',
        'obsm': 'obs'
    }
    for conflict_matrix_name, conflict in conflicts.items():
        for conflict_name, conflict_data in conflict.items():
            singledim_name = mutlidim_to_singledim[conflict_matrix_name]
            singledim_index = getattr(mod, singledim_name).index
            getattr(mod, conflict_matrix_name)[conflict_name] = conflict_data.reindex(singledim_index)
    concatenated_data.update()
    return concatenated_data

def concatenate_modalities(n_processes: int, sample_ids: tuple[str], modalities: dict[str, Iterable[anndata.AnnData]],
                           other_axis_mode: str) -> mu.MuData:
    """
    Join the modalities together into a single multimodal sample.
    """
    logger.info('Concatenating samples.')
    concat_modes = {
        "move": None,
    }
    other_axis_mode_to_apply = concat_modes.get(other_axis_mode, other_axis_mode)
    new_mods = {mod_name: anndata.concat(modes,
                                         join='outer',
                                         merge=other_axis_mode_to_apply)
                for mod_name, modes in modalities.items()}
    concatenated_data = mu.MuData(new_mods)
    logger.info('Concatenated data shape: %s', concatenated_data.shape)
    if other_axis_mode == "move":
        for mod_name, modes in modalities.items():
            conflicts, new_matrices = split_conflicts_modalities(n_processes, sample_ids, modes)
            concatenated_data = set_conflicts(concatenated_data, mod_name, conflicts)
            concatenated_data = set_matrices(concatenated_data, mod_name, new_matrices)
    logger.info("Concatenation successful.")
    return concatenated_data


def main() -> None:
    # Read in sample names and sample .h5mu files
    sample_ids: tuple[str] = tuple(i.strip() for i in par["sample_names"])
    samples: list[mu.MuData] = [mu.read(path.strip()) for path in par["input"]]

    if len(sample_ids) != len(samples):
        raise ValueError("The number of sample names must match the number of sample files.")

    if len(set(par["sample_names"])) != len(par["sample_names"]):
        raise ValueError("The sample names should be unique.")

    logger.info("\nConcatenating data for:\n\t%s\nFrom paths:\n\t%s",
                "\n\t".join(sample_ids),
                "\n\t".join(par["input"]))

    add_sample_names(sample_ids, samples, par["obs_sample_name"])
    make_observation_keys_unique(sample_ids, samples)

    mods = group_modalities(samples)
    n_processes = meta["cpus"] if meta["cpus"] else 1
    concatenated_samples = concatenate_modalities(n_processes,
                                                  sample_ids,
                                                  mods,
                                                  par["other_axis_mode"])
    logger.info("Writing out data to '%s' with compression '%s'.",
                par["output"], par["compression"])
    concatenated_samples.write(par["output"], compression=par["compression"])


if __name__ == "__main__":
    main()
