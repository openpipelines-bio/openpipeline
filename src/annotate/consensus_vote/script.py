import sys
import mudata as mu
import numpy as np
import pandas as pd

## VIASH START
par = {
    "input": "test_with_probabilities.h5mu",
    "modality": "rna",
    "input_obs_predictions": ["scanvi_pred", "celltypist_pred", "singler_pred"],
    "input_obs_probabilities": ["scanvi_prob", "celltypist_prob", "singler_prob"],
    "weights": None,
    "tie_label": None,
    "output": "consensus_test_output.h5mu",
    "output_obs_predictions": "consensus_pred",
    "output_obs_score": "consensus_score",
    "output_compression": "gzip",
}
meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from compress_h5mu import write_h5ad_to_h5mu_with_compression

logger = setup_logger()


def main():
    prediction_cols = par["input_obs_predictions"]
    prob_cols = par["input_obs_probabilities"]
    weights = par["weights"]

    if weights and len(weights) != len(prediction_cols):
        raise ValueError(
            f"--weights must have the same length as --input_obs_predictions. "
            f"Got {len(weights)} weights for {len(prediction_cols)} prediction columns."
        )
    if prob_cols and len(prob_cols) != len(prediction_cols):
        raise ValueError(
            f"--input_obs_probabilities must have the same length as --input_obs_predictions. "
            f"Got {len(prob_cols)} probability columns for {len(prediction_cols)} prediction columns."
        )

    logger.info("Reading input data.")
    adata = mu.read_h5ad(par["input"], mod=par["modality"])

    cols_to_check = [prediction_cols]
    if prob_cols:
        cols_to_check.append(prob_cols)
    for cols in cols_to_check:
        for col in cols:
            if col not in adata.obs.columns:
                raise ValueError(f"Column '{col}' not found in .obs.")

    # Each method is treated equally by default, unless user specific weights are provided
    n_methods = len(prediction_cols)
    logger.info("Initializing weights to matrix of ones")
    weights_arr = np.ones(n_methods, dtype=np.float32)
    if weights:
        logger.info("Applying user-provided weights.")
        weights_arr = np.array(weights, dtype=np.float32)
    logger.info("Normalizing weights")
    weights_arr = weights_arr / weights_arr.sum()

    # Apply the weights to the probabilities in the data
    weights = pd.DataFrame(
        [weights_arr] * adata.n_obs, index=adata.obs.index, columns=prediction_cols
    )
    if prob_cols:
        logger.info("Scaling the weights with the probabilities from each method")
        weights = weights * adata.obs[prob_cols].astype(np.float32).to_numpy()
        assert pd.notna(weights).all(axis=None)

    logger.info("Computing weighted majority vote.")
    pred_df = adata.obs[prediction_cols].astype(str)

    # For each cell and each method (index), get the label and the weight
    incidences_weights = pd.DataFrame(
        {"label": pred_df.stack(), "weights": weights.stack()}
    )
    # Move the label to the index, there might be duplicate indices now
    incidences_weights = incidences_weights.set_index("label", append=True).rename_axis(
        ["cell_id", "method", "label"]
    )
    # Sum the weights per label, from this the labels with the largest weights need to be selected
    summed_weights = incidences_weights.groupby(level=["cell_id", "label"]).sum()
    # Find the weight that is the largest per group
    max_weight_per_group = summed_weights.groupby(level="cell_id").transform("max")
    # Use the value to look-up the corresponding IDs and labels
    max_weights_mask = summed_weights["weights"] == max_weight_per_group["weights"]
    entries_for_max_weights = summed_weights[max_weights_mask].reset_index(
        level="label"
    )
    # Find the cases where there is a tie
    is_duplicated = max_weights_mask.groupby(level="cell_id").sum() > 1
    # For the ties, overwrite the label. If a cell is in the frame more than once it is because of a tie.
    entries_for_max_weights.loc[is_duplicated, ["label"]] = par["tie_label"]
    # Now its safe to just take the first index in case of duplicates, since the label and the score is the same.
    entries_for_max_weights = entries_for_max_weights[
        ~entries_for_max_weights.index.duplicated()
    ]
    # Normalize the weights
    normalized_scores = (
        entries_for_max_weights["weights"]
        / incidences_weights["weights"].groupby(level="cell_id").sum()
    )
    # Handle devision by 0
    normalized_scores = normalized_scores.replace([np.inf, -np.inf], 0.0)
    logger.info("Moving the output to the anndata.")
    adata.obs[par["output_obs_predictions"]] = entries_for_max_weights["label"].astype(
        "category"
    )
    adata.obs[par["output_obs_score"]] = normalized_scores

    logger.info("Writing output data...")
    write_h5ad_to_h5mu_with_compression(
        par["output"], par["input"], par["modality"], adata, par["output_compression"]
    )


if __name__ == "__main__":
    main()
