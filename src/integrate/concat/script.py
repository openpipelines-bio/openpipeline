from __future__ import annotations
import logging
import anndata
import muon as mu
from sys import stdout
import pandas as pd
import numpy as np

### VIASH START
par = {
    "input": ["resources_test/concat/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix_subset.h5mu",
              "resources_test/concat/human_brain_3k_filtered_feature_bc_matrix_subset.h5mu"],
    "output": "foo.h5mu",
    "sample_names": ["mouse", "human"],
    "compression": "gzip",
    "other_axis_mode": "concat"
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


def concatenate_modalities(modalities: dict[str, anndata.AnnData],
                           other_axis_mode: str) -> mu.MuData:
    """
    Join the modalities together into a single multimodal sample.
    """
    logger.info('Concatenating samples.')
    if other_axis_mode == "concat":
        other_axis_mode = concat_columns
    new_mods = {mod_name: anndata.concat(modes,
                                         join='outer',
                                         merge=other_axis_mode)
                for mod_name, modes in modalities.items()}
    concatenated_data = mu.MuData(new_mods)
    logger.info("Concatenation succesfull.")
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
