from __future__ import annotations
import logging
import anndata
import muon as mu
from sys import stdout
import pandas as pd
import numpy as np
from collections.abc import Iterable

### VIASH START
par = {
    "input": ["resources_test/concat/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix_subset.h5mu",
              "resources_test/concat/human_brain_3k_filtered_feature_bc_matrix_subset.h5mu"],
    "output": "foo.h5mu",
    "sample_names": ["mouse", "human"],
    "compression": "gzip",
    "other_axis_mode": "move"
}
### VIASH END

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)


def add_sample_names(sample_ids: tuple[str], samples: list[mu.MuData]) -> None:
    """
    Add sample names to the observations for each sample.
    Additionally, set the .batch attribute to each MuData object to store
    the sample names
    """
    for (sample_id, sample) in zip(sample_ids, samples):
        if "batch" in sample.obs_keys():
            samples.obs = sample.obs.drop("batch", axis=1)
        for _, modality in sample.mod.items():
            modality.obs["batch"] = sample_id
        sample.batch = sample_id
        for mod in sample.mod.values():
            mod.batch = sample_id
        sample.update()


def make_observation_keys_unique(samples: list[mu.MuData]) -> None:
    """
    Make the observation keys unique across all samples. At input,
    the observation keys are unique within a sample. By adding the sample name
    (unique for a sample) to each observation key, the observation key is made
    unique across all samples as well.
    """
    logger.info('Making observation keys unique across all samples.')
    for sample in samples:
        sample.obs.index = f"{sample.batch}_" + sample.obs.index
        make_observation_keys_unique_per_mod(sample)


def make_observation_keys_unique_per_mod(sample: list[anndata.AnnData]) -> None:
    """
    Updating MuData.obs_names is not allowed (it is read-only).
    So the observation keys for each modality has to be updated manually.
    """
    for _, mod in sample.mod.items():
        mod.obs_names = f"{sample.batch}_" + mod.obs_names


def group_modalities(samples: list[anndata.AnnData]) -> dict[str, anndata.AnnData]:
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


def concat_columns(vars_list: list[pd.DataFrame]) -> pd.DataFrame:
    """
    Combine dataframes by joining matching columns into a comma-separated list
    containing unique, non-na values.
    """
    column_names = set(column_name for var in vars_list for column_name in var)
    logger.debug('Trying to concatenate columns: %s.', ",".join(column_names))
    if not column_names:
        return pd.DataFrame()
    df = pd.concat(vars_list, axis=1)
    logger.debug('Concatenated dataframe. Aggregating column.')
    result = df.groupby(df.columns, axis=1).agg(lambda x: x.apply(lambda y: ','.join(y.dropna().unique().astype('str')), axis=1))
    logger.info('Aggregation completed.')
    dtypes = {}
    for col_name in column_names:
        for var in vars_list:
            col = var.get(col_name, None)
            if col is not None and col_name not in dtypes:
                dtypes[col_name] = col.dtype
    result.replace('', np.nan, inplace=True)
    result = concat_result_cast_dtype(result, orignal_dtypes=dtypes)
    logger.debug('Finished concatenating. Result is:\n%s', result)
    return result


def concat_result_cast_dtype(result: pd.DataFrame,
                             orignal_dtypes: dict[str, pd.core.dtypes.dtypes.Dtype]) -> pd.DataFrame:
    logger.debug('Trying to cast to "category" or keep original datatype.')
    for col_name, orig_dtype in orignal_dtypes.items():
        try:
            result = result.astype({col_name: "category"}, copy=True)
        except (ValueError, TypeError):
            try:
                result = result.astype({col_name: orig_dtype}, copy=True)
            except (ValueError, TypeError):
                logger.warning("Could not keep datatype for column %s", col_name)
    return result

def any_row_contains_duplicate_values(frame: pd.DataFrame) -> bool:
    """
    Check if any row contains duplicate values, that are not NA.
    """
    number_of_unique = frame.nunique(axis=1, dropna=True)
    non_na_counts = frame.count(axis=1)
    is_duplicated = (number_of_unique - non_na_counts) != 0
    return is_duplicated.any()

def split_conflicts_matrices(matrices: Iterable[pd.DataFrame], sample_ids: Iterable[str]) \
    -> tuple[dict[str, pd.DataFrame], pd.DataFrame | None]:
    """
    Merge matrices by combining columns that have the same name.
    Columns that contain conflicting values (e.i. the columns have different values),
    are not merged, but instead moved to a new dataframe.
    """
    column_names = set(column_name for var in matrices for column_name in var)
    logger.debug('Trying to concatenate columns: %s.', ",".join(column_names))
    if not column_names:
        return {}, None
    conflicts = {}
    concatenated_matrix = pd.DataFrame()
    for column_name in column_names:
        columns = [var[column_name] for var in matrices if column_name in var]
        assert columns, "Some columns should have been found."
        concatenated_columns = pd.concat(columns, axis=1)
        if any_row_contains_duplicate_values(concatenated_columns):
            concatenated_columns.columns = sample_ids
            conflicts[f'conflict_{column_name}'] = concatenated_columns
        else:
            unique_values = concatenated_columns.fillna(method='bfill', axis=1).iloc[:, 0]
            concatenated_matrix = concatenated_matrix.assign(**{column_name: unique_values})
    return conflicts, concatenated_matrix

def split_conflicts_modalities(modalities: Iterable[anndata.AnnData]) \
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
            sample_ids = [modality.batch for modality in modalities]
            conflicts, concatenated_matrix = split_conflicts_matrices(matrices, sample_ids)
            conflicts_result[f"{matrix_name}m"] = conflicts
            concatenated_result[matrix_name] = concatenated_matrix
        return conflicts_result, concatenated_result

def set_matrices(concatenated_data: mu.MuData, 
                 mod_name: str,
                 new_matrices: dict[str, pd.DataFrame | None]) -> mu.MuData:
    mod = concatenated_data.mod[mod_name]
    for matrix_name, data in new_matrices.items():
        if data is None:
            data = pd.DataFrame(index=getattr(mod, matrix_name).index)
        setattr(mod, matrix_name, data)
    return concatenated_data

def set_conflicts(concatenated_data: mu.MuData,
                  mod_name: str,
                  conflicts: dict[str, dict[str, pd.DataFrame]]) -> mu.MuData:
    mod = concatenated_data.mod[mod_name]
    for conflict_matrix_name, conflict in conflicts.items():
        for conflict_name, conflict_data in conflict.items():
            getattr(mod, conflict_matrix_name)[conflict_name] = conflict_data.sort_index()
    return concatenated_data

def concatenate_modalities(modalities: dict[str, anndata.AnnData],
                           other_axis_mode: str) -> mu.MuData:
    """
    Join the modalities together into a single multimodal sample.
    """
    logger.info('Concatenating samples.')
    concat_modes = {
        "concat": concat_columns,
        "move": "same",
    }
    other_axis_mode_to_apply = concat_modes.get(other_axis_mode, other_axis_mode)
    new_mods = {mod_name: anndata.concat(modes,
                                         join='outer',
                                         merge=other_axis_mode_to_apply)
                for mod_name, modes in modalities.items()}
    concatenated_data = mu.MuData(new_mods)
    if other_axis_mode == "move":
        for mod_name, modes in modalities.items():
            conflicts, new_matrices = split_conflicts_modalities(modes)
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

    add_sample_names(sample_ids, samples)
    make_observation_keys_unique(samples)

    mods = group_modalities(samples)
    concatenated_samples = concatenate_modalities(mods, par["other_axis_mode"])
    logger.info("Writing out data to '%s' with compression '%s'.",
                par["output"], par["compression"])
    concatenated_samples.write(par["output"], compression=par["compression"])


if __name__ == "__main__":
    main()
