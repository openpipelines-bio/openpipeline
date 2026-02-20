import sys
import h5py
import mudata
from anndata.io import read_elem, write_elem
from anndata import AnnData, settings
from scipy.sparse import csr_array
import numpy as np
from contextlib import contextmanager, closing, nullcontext
from shutil import copytree, copyfile
from functools import partial
import zarr

settings.zarr_write_format = 3
# Avoid FutureWarning for mudata 0.4
mudata.set_options(pull_on_update=True)

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "output": "foo.h5mu",
    "modality": "rna",
    "layer": None,
    "top_n_vars": [10, 20, 50],
    "var_qc_metrics": None,
    "output_var_obs_mean": "obs_mean",
    "output_var_total_counts_obs": "total_counts",
    "output_var_num_nonzero_obs": "num_nonzero_obs",
    "output_var_pct_dropout": "pct_dropout",
    "output_obs_total_counts_vars": "total_counts",
    "output_obs_num_nonzero_vars": "num_nonzero_vars",
    "output_compression": None,
}
meta = {"resources_dir": "./src/utils", "name": "calculate_qc_metrics"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()


@contextmanager
def mudata_opener(file_loc, mode=None):
    open_mudata = None
    input_is_zarr = False
    try:
        open_mudata = zarr.open(file_loc, zarr_format=3, mode=mode)
        input_is_zarr = True
        yield open_mudata, input_is_zarr
    except (zarr.errors.GroupNotFoundError, NotADirectoryError):
        try:
            open_mudata = h5py.File(file_loc, mode=mode)
            yield open_mudata, input_is_zarr
        except (FileNotFoundError, IsADirectoryError, KeyError) as e:
            e.add_note(f"Could not open file {file_loc}.")
            raise e
        finally:
            try:
                if open_mudata:
                    open_mudata.close()
                    del open_mudata
            except (AttributeError, UnboundLocalError):
                pass


def calculate_var_statistics(layer):
    logger.info("Calculating statistics to store in .var")
    var_columns_to_add = {}
    if par["output_var_obs_mean"]:
        logger.info(
            "(var) Calculating mean per observation, which will be stored at %s.",
            par["output_var_obs_mean"],
        )
        obs_mean = np.ravel(mean_csr_array(layer, axis=0))
        var_columns_to_add[par["output_var_obs_mean"]] = obs_mean
    if par["output_var_total_counts_obs"]:
        logger.info(
            "(var) Calculating total counts for each observation, to be stored at %s.",
            par["output_var_total_counts_obs"],
        )
        total_counts_obs = np.ravel(layer.sum(axis=0))
        var_columns_to_add[par["output_var_total_counts_obs"]] = total_counts_obs

    # This is the same as the old .nnz(axis=0), but this new implementation only works for csr_arrays!
    # See https://github.com/scipy/scipy/issues/19405#issuecomment-1773553180
    num_nonzero_obs = np.bincount(layer.indices, minlength=layer.shape[1])
    if par["output_var_num_nonzero_obs"]:
        logger.info(
            "(var) Retreiving the number of non-zero elements for each row for result column %s",
            par["output_var_num_nonzero_obs"],
        )
        var_columns_to_add[par["output_var_num_nonzero_obs"]] = num_nonzero_obs
    if par["output_var_pct_dropout"]:
        logger.info(
            "(var) Fetching for each feature the percentage of observations missing that feature (column %s)",
            par["output_var_pct_dropout"],
        )
        var_columns_to_add[par["output_var_pct_dropout"]] = (
            1 - num_nonzero_obs / layer.shape[0]
        ) * 100
    logger.info("Calculating .var statistics finished.")
    return var_columns_to_add


def calculate_obs_statistics(layer, var):
    logger.info("Calculating statistics to store in .obs")
    obs_columns_to_add = {}
    total_counts_var = np.ravel(layer.sum(axis=1))

    if par["output_obs_num_nonzero_vars"]:
        logger.info(
            "(obs) Retreiving the number of non-zero elements for each feature to be stored in column %s",
            par["output_var_num_nonzero_obs"],
        )
        # This is the same as the old .nnz(axis=1), but this new implementation only works for csr_arrays!
        # See https://github.com/scipy/scipy/issues/19405#issuecomment-1773553180
        num_nonzero_vars = np.diff(layer.indptr)
        obs_columns_to_add[par["output_obs_num_nonzero_vars"]] = num_nonzero_vars

    if par["output_obs_total_counts_vars"]:
        logger.info(
            "(obs) Calculating total counts for each feature, to be stored at %s.",
            par["output_obs_total_counts_vars"],
        )
        obs_columns_to_add[par["output_obs_total_counts_vars"]] = total_counts_var

    top_metrics = {}
    if par["top_n_vars"]:
        logger.info(
            "(obs) Calculating the cumulative proportions to the %s most expressed vars.",
            ", ".join([f"{top_n_var}th" for top_n_var in par["top_n_vars"]]),
        )
        par["top_n_vars"] = sorted(par["top_n_vars"])
        distributions = get_top_from_csr_matrix(layer, par["top_n_vars"])
        top_metrics = {
            distribution_size: distribution * 100
            for distribution_size, distribution in zip(
                par["top_n_vars"], distributions.T
            )
        }
        obs_columns_to_add |= {
            f"pct_of_counts_in_top_{n_top}_vars": col
            for n_top, col in top_metrics.items()
        }
    for qc_metric in par.get("var_qc_metrics", []) or []:
        logger.info(
            "(obs) Retreiving the proportion of total 'True' values in column %s",
            qc_metric,
        )
        if qc_metric not in var:
            raise ValueError(
                f"Value for --var_qc_metrics, '{qc_metric}' "
                f"not found in .var for modality {par['modality']}"
            )
        qc_column = var[qc_metric]
        if qc_column.isna().any():
            if par["var_qc_metrics_fill_na_value"] is None:
                raise ValueError(
                    f"The .var column '{qc_metric}', selected by '--var_qc_metrics', contains NA values. "
                    "It is ambiguous whether or not to include these values in the static calulation. "
                    "You can explicitly map the NA values to 'False' or 'True using '--var_qc_metrics_fill_na_value'"
                )
            else:
                qc_column = qc_column.fillna(
                    par["var_qc_metrics_fill_na_value"], inplace=False
                )
        qc_column = qc_column.to_list()
        if set(np.unique(qc_column)) - {True, False}:
            raise ValueError(
                f"Column {qc_metric} in .var for modality {par['modality']} "
                f"must only contain boolean values"
            )
        total_counts_qc_metric = np.ravel(layer[:, qc_column].sum(axis=1))
        obs_columns_to_add |= {
            f"total_counts_{qc_metric}": total_counts_qc_metric,
            f"pct_{qc_metric}": total_counts_qc_metric / total_counts_var * 100,
        }
    logger.info("Finised calculating obs statistics")
    return obs_columns_to_add


def cast_layer_dtype(layer):
    # from the np.sum documentation:
    # Especially when summing a large number of lower precision floating point numbers,
    # such as float32, numerical errors can become significant. In such cases it can
    # be advisable to use dtype="float64" to use a higher precision for the output.

    # However, the 'dtype' from SciPy's implementation of sum cannot be used for this
    # as it does not use an internal accumulator but matrix multiplication for calculating
    # the sum. So here we cast explicitly to a higher precision before doing the calcualtions.
    # The downside is that this is inefficient.
    # See https://github.com/scipy/scipy/issues/23768#issuecomment-3909317463
    original_dtype = layer.dtype
    target_dtype = original_dtype
    if np.issubdtype(original_dtype, np.floating) and np.can_cast(
        original_dtype, np.float64, casting="safe"
    ):
        # use promote_types in orde to make suresure not to cast np.float128
        # or anything else to a lower precision dtype
        target_dtype = np.promote_types(np.float64, original_dtype)
        logger.info(
            "Using target dtype %s for layer. Casting may be required for higher precision.",
            target_dtype,
        )
    result = csr_array(layer, dtype=target_dtype, copy=False)
    logger.info("Constructed CSR of shape %s and dtype %s", result.shape, result.dtype)
    return result


def mean_csr_array(input_csr_array, axis):
    # TODO: replace this function with the native SciPy version
    # when anndata supports SciPy >= 1.17.0
    # We use this version to avoid creating duplicates of the data in memory
    # when the `astype` function is called.
    # See https://github.com/scipy/scipy/pull/23797
    sparse_array_dtype = input_csr_array.dtype
    integral = np.issubdtype(sparse_array_dtype, np.integer) or np.issubdtype(
        sparse_array_dtype, np.bool_
    )

    # intermediate dtype for summation
    inter_dtype = np.float64 if integral else sparse_array_dtype
    inter_cast = input_csr_array.astype(inter_dtype, copy=False)
    divided = inter_cast.data * (1.0 / input_csr_array.shape[axis])
    divided_csr = csr_array(
        (divided, input_csr_array.indices, input_csr_array.indptr),
        shape=input_csr_array.shape,
    )
    return divided_csr.sum(axis=axis, dtype=inter_dtype)


def get_top_from_csr_matrix(array, top_n_genes):
    # csr matrices stores a 3D matrix in a format such that data for individual cells
    # are stored in 1 array. Another array (indptr) here stores the ranges of indices
    # to select from the data-array (.e.g. data[indptr[0]:indptr[1]] for row 0) for each row.
    # Another array 'indices' maps each element of data to a column
    # (data and indices arrays have the same length)
    top_n_genes = np.array(top_n_genes).astype(np.int64)
    assert np.all(top_n_genes[:-1] <= top_n_genes[1:]), "top_n_genes must be sorted"
    row_indices, data = array.indptr, array.data
    number_of_rows, max_genes_to_parse = row_indices.size - 1, top_n_genes[-1]
    top_data = np.zeros((number_of_rows, max_genes_to_parse), dtype=data.dtype)
    # Loop over each row to create a dense matrix without the 0 counts,
    # but not for the whole matrix, only store the genes up until
    # the largest number of top n genes.
    for row_number in range(number_of_rows):
        row_start_index, row_end_index = (
            row_indices[row_number],
            row_indices[row_number + 1],
        )
        row_data = data[row_start_index:row_end_index]  # all non-zero counts for an row
        try:
            # There are less genes with counts in the row than the
            # maximum number of genes we would like to select
            # all these genes are in the top genes, just store them
            top_data[row_number, : row_end_index - row_start_index] = row_data
        except ValueError:
            # Store the counts for the top genes
            top_data[row_number, :] = np.partition(row_data, -max_genes_to_parse)[
                -max_genes_to_parse:
            ]

    # Partition works from smallest to largest, but we want largest
    # so do smallest to largest first (but with reversed indices)
    top_data = np.partition(top_data, max_genes_to_parse - top_n_genes)
    # And then switch the order around
    top_data = np.flip(top_data, axis=1)

    cumulative = top_data.cumsum(axis=1, dtype=np.float64)[:, top_n_genes - 1]
    return cumulative / np.expand_dims(array.sum(axis=1), 1)


def main():
    logger.info("Started %s.", meta["name"])
    mod_element_loc = f"mod/{par['modality']}"
    layer_element_name = (
        f"{mod_element_loc}/X"
        if not par["layer"]
        else f"{mod_element_loc}/layers/{par['layer']}"
    )
    # In order to match the format (zarr or h5) of the input to the output
    modality_anndatas = {}
    with mudata_opener(par["input"], mode="r") as (open_mudata, input_is_zarr):
        logger.info(
            "Openened %s in %s format.", par["input"], "zarr" if input_is_zarr else "h5"
        )
        mods = list(open_mudata["mod"].keys())
        logger.info("Found modalities: %s", ", ".join(mods))
        # Create AnnData object for only the metadata (var and obs)
        # We need all of the modalities for this in order to be able to adjust
        # the 'global' (across all modalities) dataframes in the MuData object
        # when writing back the output.
        for mod in mods:
            logger.info("Reading metadata frames for modality %s", mod)
            var = read_elem(open_mudata[f"/mod/{mod}/var"])
            logger.info(".var shape for %s is %s", mod, var.shape)
            obs = read_elem(open_mudata[f"/mod/{mod}/obs"])
            logger.info(".obs shape for %s is %s", mod, obs.shape)
            modality_anndatas[mod] = AnnData(var=var, obs=obs)

        logger.info("Reading layer %s", "X" if not par["layer"] else par["layer"])
        layer = read_elem(open_mudata[layer_element_name])
        logger.info("Found layer with shape %s and dtype %s", layer.shape, layer.dtype)

    layer = cast_layer_dtype(layer)
    logger.info("Eliminating explicit zeros from sparse layer.")
    layer.eliminate_zeros()

    var_columns_to_add = calculate_var_statistics(layer)
    modality_anndatas[par["modality"]].var = modality_anndatas[
        par["modality"]
    ].var.assign(**var_columns_to_add)

    # obs statistics
    obs_columns_to_add = calculate_obs_statistics(
        layer, modality_anndatas[par["modality"]].var
    )
    modality_anndatas[par["modality"]].obs = modality_anndatas[
        par["modality"]
    ].obs.assign(**obs_columns_to_add)

    # Use MuData internals to construct global metadata dataframes
    logger.info("Constructing global (across all modalities) obs and var dataframes")
    mudata_skeleton = mudata.MuData(modality_anndatas)
    logger.info(
        "Global var and obs dataframes had shape %s and %s respectively",
        mudata_skeleton.var.shape,
        mudata_skeleton.obs.shape,
    )

    logger.info("Writing to %s", par["output"])
    try:
        copytree(par["input"], par["output"], symlinks=True)
    except NotADirectoryError:
        copyfile(par["input"], par["output"], follow_symlinks=True)
        logger.info("Copied input file %s to %s", par["input"], par["output"])
    else:
        logger.info("Copied input directory %s", par["input"], par["output"])
    logger.info("Using a %s writer", "zarr" if input_is_zarr else "H5")
    write_opener = partial(zarr.open, zarr_format=3) if input_is_zarr else h5py.File
    context = (
        nullcontext if input_is_zarr else closing
    )  # zarr format does not need to be closed
    logger.info("Overwriting slots.")
    with context(write_opener(par["output"], mode="a")) as open_output:
        write_elem(
            open_output[mod_element_loc], "obs", modality_anndatas[par["modality"]].obs
        )
        write_elem(
            open_output[mod_element_loc], "var", modality_anndatas[par["modality"]].var
        )
        write_elem(open_output, "var", mudata_skeleton.var)
        write_elem(open_output, "obs", mudata_skeleton.obs)
    logger.info("Finished!")


if __name__ == "__main__":
    main()
