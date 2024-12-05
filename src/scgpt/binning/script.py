import sys
import mudata as mu
import numpy as np
from scipy.sparse import csr_matrix
import warnings

## VIASH START
par = {
    "input": "resources_test/scgpt/test_resources/Kim2020_Lung_subset_genes_cross_checked.h5mu",
    "output": "resources_test/scgpt/test_resources/Kim2020_Lung_subset_binned.h5mu",
    "modality": "rna",
    "input_layer": None,
    "output_obsm_binned_counts": "binned_counts",
    "n_input_bins": 51,
    "output_compression": None,
    "var_input": "id_in_vocab",
    "seed": 0
}
meta = {
    "resources_dir": "src/utils"
}
## VIASH END

if par["seed"]:
    np.random.seed(par["seed"])

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from subset_vars import subset_vars
logger = setup_logger()

logger.info("Reading in data")
# Read in data
mdata = mu.read(par["input"])
input_adata = mdata.mod[par["modality"]]
adata = input_adata.copy()

logger.info("Subsetting data based on highly variable gene and/or cross-checked genes")
adata = subset_vars(adata, par["var_input"])

logger.info("Converting the input layer into a CSR matrix")
if not par['input_layer'] or par["input_layer"] == "X":
    layer_data = adata.X
else:
    layer_data = adata.layers[par['input_layer']]
layer_data = csr_matrix(layer_data)

if layer_data.min() < 0:
    raise ValueError(
        f"Assuming non-negative data, but got min value {layer_data.min()}."
    )

n_bins = par["n_input_bins"]  # NOTE: the first bin is always a spectial for zero
logger.info(f"Binning data into {par['n_input_bins']} bins.")


def _digitize(x: np.ndarray, bins: np.ndarray) -> np.ndarray:
    assert x.ndim == 1 and bins.ndim == 1

    left_digits = np.digitize(x, bins)
    right_difits = np.digitize(x, bins, right=True)

    rands = np.random.rand(len(x))  # uniform random numbers

    digits = rands * (right_difits - left_digits) + left_digits
    digits = np.ceil(digits)
    smallest_dtype = np.min_scalar_type(digits.max().astype(np.uint)) # Already checked for non-negative values
    digits = digits.astype(smallest_dtype)

    return digits


with warnings.catch_warnings():
    # Make sure warnings are displayed once.
    warnings.simplefilter("once")
    # layer_data.indptr.size is the number of rows in the sparse matrix 
    binned_rows = []
    bin_edges = []
    logger.info("Establishing bin edges and digitizing of non-zero values into bins for each row of the count matrix")
    for row_number in range(layer_data.indptr.size-1):
        row_start_index, row_end_index = layer_data.indptr[row_number], layer_data.indptr[row_number+1]
        # These are all non-zero counts in the row
        non_zero_row = layer_data.data[row_start_index:row_end_index] 
        if len(non_zero_row) == 0:
            logger.warning(
                "The input data contains all zero rows. Please make sure "
                "this is expected. You can use the `filter_cell_by_counts` "
                "arg to filter out all zero rows."
            )

            # Add binned_rows and bin_edges as all 0
            # np.stack will upcast the dtype later
            binned_rows.append(np.zeros_like(non_zero_row, dtype=np.int8))
            bin_edges.append(np.array([0] * n_bins))
            continue

        # Binning of non-zero values
        bins = np.quantile(non_zero_row, np.linspace(0, 1, n_bins - 1))
        non_zero_digits = _digitize(non_zero_row, bins)
        assert non_zero_digits.min() >= 1
        assert non_zero_digits.max() <= n_bins - 1
        binned_rows.append(non_zero_digits)

        bin_edges.append(np.concatenate([[0], bins]))

# Create new CSR matrix
logger.info("Creating a new CSR matrix of the binned count values")
binned_counts = csr_matrix((np.concatenate(binned_rows, casting="same_kind"), 
                          layer_data.indices, layer_data.indptr), shape=layer_data.shape)

# Set binned values and bin edges layers to adata object
input_adata.obsm[par["output_obsm_binned_counts"]] = binned_counts
input_adata.obsm["bin_edges"] = np.stack(bin_edges)

# Write mudata output 
logger.info("Writing output data")
mdata. write_h5mu(par["output"], compression=par["output_compression"])
