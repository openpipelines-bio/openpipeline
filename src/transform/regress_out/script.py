import scanpy as sc
import mudata as mu
import anndata as ad
import multiprocessing
import sys
import numpy as np
from scipy.sparse import csr_matrix

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "obs_keys": [],
    "var_input": None,
}
meta = {"name": "lognorm"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

logger.info("Reading input mudata")
mdata = mu.read_h5mu(par["input"])
mdata.var_names_make_unique()

if par["obs_keys"] is not None and len(par["obs_keys"]) > 0:
    mod = par["modality"]
    data = mdata.mod[mod]

    # Copy required data from input data to new AnnData object to allow providing input and output layers
    input_layer = data.X if not par["input_layer"] else data.layers[par["input_layer"]]

    mask_var = None
    if par["var_input"]:
        mask_var = data.var[par["var_input"]].to_numpy(dtype=bool)
        logger.info(
            "Regressing out on %i genes selected by .var column %s",
            mask_var.sum(),
            par["var_input"],
        )

    obs = data.obs.loc[:, par["obs_keys"]]
    X = input_layer.copy() if mask_var is None else input_layer[:, mask_var]
    sc_data = ad.AnnData(X=X, obs=obs)

    logger.info("Regress out variables on modality %s", mod)
    sc.pp.regress_out(
        sc_data, keys=par["obs_keys"], n_jobs=multiprocessing.cpu_count() - 1
    )

    if mask_var is None:
        regressed = sc_data.X
    else:
        # Store only the selected genes in a sparse matrix; non-selected genes are 0.
        regressed_subset = csr_matrix(sc_data.X)
        del sc_data, X
        selected_indices = np.flatnonzero(mask_var).astype(
            regressed_subset.indices.dtype
        )
        regressed = csr_matrix(
            (
                regressed_subset.data,
                selected_indices[regressed_subset.indices],
                regressed_subset.indptr,
            ),
            shape=(data.n_obs, data.n_vars),
            copy=False,
        )
        del regressed_subset

    # Copy regressed data back to original input data
    if par["output_layer"]:
        data.layers[par["output_layer"]] = regressed
    else:
        data.X = regressed

logger.info("Writing to file")
mdata.write_h5mu(filename=par["output"], compression=par["output_compression"])
