import numpy as np
import mudata as mu
import pandas as pd
import sys
import scanpy as sc
import scipy.sparse as sp

## VIASH START
par = {
    "input": "resources_test/annotation_test_data/TS_Blood_filtered.h5mu",
    "modality": "rna",
    "input_layer": None,
    "obs_label": "cell_type",
    "obs_groups": ["treatment", "donor_id", "disease"],
    "obs_cell_count": "n_cells",
    "random_state": 0,
    "output": "test.h5mu",
    "output_compression": "gzip",
}
meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()


def is_normalized(layer):
    exp_layer = np.expm1(layer)  # Inverse of log1p

    if sp.issparse(layer):
        row_sums = np.array(layer.sum(axis=1)).flatten()
        exp_row_sums = np.array(exp_layer.sum(axis=1)).flatten()
    else:
        row_sums = layer.sum(axis=1)
        exp_row_sums = exp_layer.sum(axis=1)

    is_normalized = np.allclose(row_sums, 1)
    is_log1p_normalized = np.isfinite(exp_row_sums).all() and np.allclose(
        exp_row_sums, exp_row_sums[0]
    )

    return is_normalized or is_log1p_normalized


def count_obs(adata, pb_adata, obs_cols):
    counts = []
    for i in range(pb_adata.n_obs):
        values = pb_adata.obs.iloc[i][obs_cols]

        mask = pd.Series([True] * adata.n_obs)
        for col in obs_cols:
            mask &= (adata.obs[col] == values[col]).values

        count = mask.sum()
        counts.append(count)

    return counts


def main():
    # Read in data
    logger.info(f"Reading input data {par['input']}...")
    adata = mu.read_h5ad(par["input"], mod=par["modality"]).copy()

    # Make sure .X contains raw counts
    if par["input_layer"]:
        adata.X = adata.layers[par["input_layer"]]
    if is_normalized(adata.X):
        raise ValueError("Input layer must contain raw counts.")

    # Sanitize pseudobulk aggregation fields
    pseudobulk_cols = [par["obs_label"]]
    if par["obs_groups"]:
        pseudobulk_cols += par["obs_groups"]

    for col in pseudobulk_cols:
        if col not in adata.obs.columns:
            raise ValueError(f"Required column '{col}' not found in .obs.")
        adata.obs[col] = (
            adata.obs[col]
            .astype(str)
            .str.replace(" ", "_")
            .str.replace("+", "")
            .astype("category")
        )

    # Aggregate pseudobulk data per cell group
    logger.info("Creating pseudobulk samples...")
    adata_pb = sc.get.aggregate(adata, by=pseudobulk_cols, func="sum", axis=0)
    adata_pb.X = adata_pb.layers["sum"]
    del adata_pb.layers["sum"]

    adata_pb.obs_names = [
        "_".join(values)
        for values in zip(*[adata_pb.obs[col] for col in pseudobulk_cols])
    ]

    # Filter pseudobulk samples based on minimum observation count
    logger.info("Filtering pseudobulk samples based on minimum observation count...")
    adata_pb.obs[par["obs_cell_count"]] = count_obs(adata, adata_pb, pseudobulk_cols)

    logger.info(
        f"Final dataset: {adata_pb.n_obs} pseudobulk samples, {adata_pb.n_vars} genes"
    )

    logger.info("Writing output data...")
    mdata_pb = mu.MuData({"rna": adata_pb})
    mdata_pb.write_h5mu(par["output"], compression=par["output_compression"])


if __name__ == "__main__":
    main()
