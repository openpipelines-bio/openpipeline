import mudata as mu
import numpy as np
from scanpy.get import _get_obs_rep
from sklearn.utils import issparse
from scipy.sparse import csr_matrix
import warnings

## VIASH START
par = {
    "input": "resources_test/scgpt/test_resources/Kim2020_Lung_subset.h5mu",
    "output": "resources_test/scgpt/test_resources/Kim2020_Lung_subset_binned_sparse.h5mu",
    "modality": "rna",
    "input_layer": "X",
    "binned_layer": "X_binned",
    "n_input_bins": 51,
    "output_compression": None
}
## VIASH END
np.random.seed(0)

# START TEMPORARY WORKAROUND setup_logger
# reason: resources aren't available when using Nextflow fusion
# from setup_logger import setup_logger
def setup_logger():
    import logging
    from sys import stdout

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    console_handler = logging.StreamHandler(stdout)
    logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
    console_handler.setFormatter(logFormatter)
    logger.addHandler(console_handler)

    return logger
# END TEMPORARY WORKAROUND setup_logger
logger = setup_logger()

logger.info("Reading in data")
# Read in data
mdata = mu.read(par["input"])
input_adata = mdata.mod[par["modality"]]
adata = input_adata.copy()


logger.info(f"Binning data of {par['input_layer']} layer into {par['n_input_bins']} bins.")
key_to_process = par["input_layer"]
if key_to_process == "X":
    key_to_process = None # Scanpy API uses arg None to use X
    
n_bins = par["n_input_bins"]  # NOTE: the first bin is always a spectial for zero

binned_rows = []
bin_edges = []

layer_data = _get_obs_rep(adata, layer=key_to_process)
layer_data = csr_matrix(layer_data)

if layer_data.min() < 0:
    raise ValueError(
        f"Assuming non-negative data, but got min value {layer_data.min()}."
    )


def _digitize(x: np.ndarray, bins: np.ndarray, side="both") -> np.ndarray:
    assert x.ndim == 1 and bins.ndim == 1

    left_digits = np.digitize(x, bins)
    if side == "one":
        return left_digits

    right_difits = np.digitize(x, bins, right=True)

    rands = np.random.rand(len(x))  # uniform random numbers

    digits = rands * (right_difits - left_digits) + left_digits
    digits = np.ceil(digits).astype(np.int64)
    return digits


with warnings.catch_warnings():
    # Make sure warnings are displayed once.
    warnings.simplefilter("once")
    # layer_data.indptr.size is the number of rows in the sparse matrix 
    for row_number in range(layer_data.indptr.size-1):
        row_start_index, row_end_index = layer_data.indptr[row_number], layer_data.indptr[row_number+1]
        # These are all non-zero counts in the row
        non_zero_row = layer_data.data[row_start_index:row_end_index] 
        if non_zero_row.max() == 0:
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
binned_layer = csr_matrix((np.concatenate(binned_rows, casting="same_kind"), 
                          layer_data.indices, layer_data.indptr), shape=layer_data.shape)

# Set binned values and bin edges layers to adata object
adata.layers[par["binned_layer"]] = binned_layer 
adata.obsm["bin_edges"] = np.stack(bin_edges)
      
# Write mudata output 
mdata.mod[par["modality"]] = adata
mdata.write(par["output"], compression=par["output_compression"]) 