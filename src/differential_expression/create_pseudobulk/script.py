import random
import numpy as np
import anndata as ad
import mudata as mu
import logging
import sys
import scipy.sparse as sp

## VIASH START
adata = mu.read_h5ad(
    "resources_test/annotation_test_data/TS_Blood_filtered.h5mu", mod="rna"
)
np.random.seed(0)
n_cells = adata.n_obs
treatment = np.random.choice(["ctrl", "stim"], size=n_cells, p=[0.5, 0.5])
disease = np.random.choice(["healthy", "diseased"], size=n_cells, p=[0.5, 0.5])
adata.obs["treatment"] = treatment
adata.obs["disease"] = disease
mdata = mu.MuData({"rna": adata})
mdata.write_h5mu("resources_test/annotation_test_data/TS_Blood_filtered_annotated.h5mu")

par = {
    "input": "resources_test/annotation_test_data/TS_Blood_filtered_annotated.h5mu",
    "modality": "rna",
    "input_layer": None,
    "obs_grouping": "cell_type",
    "obs_sample_conditions": ["treatment", "donor_id", "disease"],
    # "obs_sample_conditions": ["donor_id", "treatment", "disease"],
    "obs_pseudo_bulk_samples": "pb_sample",
    "aggregation_method": "sum",
    "min_num_cells_per_sample": 5,
    "pseudo_replicates": 1,
    "random_state": 0,
    "output": "resources_test/annotation_test_data/TS_Blood_filtered_annotated_pseudobulk.h5mu",
}
meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()


def aggregate_and_filter(
    adata,
    grouping_value,
    group_key,
    sample_key,
    obs_to_keep,
    aggregation_method="sum",
    num_cells_per_donor=30,
    pseudo_replicates=1,
):
    # subset adata to the given grouping (cell/cluster) value
    adata_cell_pop = adata[adata.obs[group_key] == grouping_value]

    pb_aggregated_data = []
    for i, pbs in enumerate(
        pb_samples := adata_cell_pop.obs[sample_key].cat.categories
    ):
        logger.info(f"Processing pseudobulk sample {i + 1} out of {len(pb_samples)}...")

        pb_adata = adata_cell_pop[adata_cell_pop.obs[sample_key] == pbs]

        # check which pseudobulk sample to keep according to the number of cells
        if pb_adata.shape[0] < num_cells_per_donor:
            logging.info(
                f"Number of cells for {pbs} is less than {num_cells_per_donor}, skipping pseudobulk aggregation..."
            )
            continue

        # Create pseudo-replicates
        logging.info(
            f"Creating {pseudo_replicates} pseudo-replicates for pseudobulk sample {pbs}..."
        )
        indices = list(pb_adata.obs_names)
        np.random.seed(par["random_state"])  # optional seed for reproducibility
        random.shuffle(indices)
        indices = np.array_split(np.array(indices), pseudo_replicates)

        # Aggregate data for each pseudo-replicate
        for i, rep_idx in enumerate(indices):
            # Aggregate gene expression values
            if aggregation_method == "sum":
                pb_adata_rep = ad.AnnData(
                    X=pb_adata[rep_idx].X.sum(axis=0), var=pb_adata[rep_idx].var[[]]
                )
            elif aggregation_method == "mean":
                pb_adata_rep = ad.AnnData(
                    X=pb_adata[rep_idx].X.mean(axis=0), var=pb_adata[rep_idx].var[[]]
                )

            # Store metadata
            pb_adata_rep.obs_names = [f"{grouping_value}_{pbs}_{str(i)}"]
            pb_adata_rep.obs[group_key] = grouping_value
            for obs in obs_to_keep:
                pb_adata_rep.obs[obs] = pb_adata[rep_idx].obs[obs][0]
            pb_adata_rep.obs["n_cells"] = len(pb_adata[rep_idx])
            pb_aggregated_data.append(pb_adata_rep)

    adata_pb = ad.concat(pb_aggregated_data)

    return adata_pb


def is_normalized(layer):
    if sp.issparse(layer):
        row_sums = np.array(layer.sum(axis=1)).flatten()
    else:
        row_sums = layer.sum(axis=1)

    return np.allclose(row_sums, 1)


def main():
    # Read in data
    logger.info(f"Reading input data {par['input']}...")
    adata = mu.read_h5ad(par["input"], mod=par["modality"]).copy()

    # Make sure .X contains raw counts
    if par["input_layer"]:
        adata.X = adata.layers[par["input_layer"]]
    if is_normalized(adata.X):
        raise ValueError("Input layer must contain raw counts.")

    # Ensure al columns required for pseudobulk aggregation exist and are categorical
    required_categorical_obs = par["obs_sample_conditions"] + [par["obs_grouping"]]
    for col in required_categorical_obs:
        if col not in adata.obs.columns:
            raise ValueError(f"Required column '{col}' not found in .obs.")
        adata.obs[col] = (
            adata.obs[col]
            .astype(str)
            .replace(" ", "_")
            .replace("+", "")
            .astype("category")
        )

    # Create pseudobulk groups
    logger.info("Creating pseudobulk groups...")
    adata.obs[par["obs_pseudo_bulk_samples"]] = [
        "_".join(values)
        for values in zip(*[adata.obs[col] for col in par["obs_sample_conditions"]])
    ]
    adata.obs[par["obs_pseudo_bulk_samples"]] = adata.obs[
        par["obs_pseudo_bulk_samples"]
    ].astype("category")
    logger.info(
        f"Pseudobulk groups created: {adata.obs[par['obs_pseudo_bulk_samples']].unique().tolist()}"
    )

    # Aggregate pseudobulk data per cell group
    pb_datasets = []
    cell_groups = adata.obs[par["obs_grouping"]].unique()

    for group in cell_groups:
        logger.info(f"Processing cell group {group}")
        pb_data = aggregate_and_filter(
            adata,
            group,
            par["obs_grouping"],
            par["obs_pseudo_bulk_samples"],
            par["obs_sample_conditions"],
            aggregation_method=par["aggregation_method"],
            num_cells_per_donor=par["min_num_cells_per_sample"],
            pseudo_replicates=par["pseudo_replicates"],
        )
        pb_datasets.append(pb_data)

    # Combine pseudobulk datasets
    logger.info("Combining pseudobulk datasets...")
    adata_pb = ad.concat(pb_datasets)
    logger.info(
        f"Final dataset: {adata_pb.n_obs} pseudobulk samples, {adata_pb.n_vars} genes"
    )
    logger.info("Writing output data...")
    mdata_pb = mu.MuData({"rna": adata_pb})
    mdata_pb.write_h5mu(par["output"], compression=par["output_compression"])


if __name__ == "__main__":
    main()
