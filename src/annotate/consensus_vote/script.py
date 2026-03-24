import sys
import mudata as mu
import numpy as np
import pandas as pd

## VIASH START
par = {
    "input": "consensus_test.h5mu",
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


def _weighted_majority(row, weights, tie_label=None):
    weights = weights.loc[row.name]
    vote_totals = {}
    for label, w in zip(row, weights):
        vote_totals[label] = vote_totals.get(label, 0.0) + w
    total = sum(vote_totals.values())
    top_score = max(vote_totals.values())
    winners = [label for label, score in vote_totals.items() if score == top_score]
    normalized_score = top_score / total if total > 0 else 0.0
    if len(winners) > 1:
        return tie_label, normalized_score
    return winners[0], normalized_score


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

    n_methods = len(prediction_cols)
    weights_arr = (
        np.array(weights, dtype=np.float32)
        if weights
        else np.ones(n_methods, dtype=np.float32)
    )
    weights_arr = weights_arr / weights_arr.sum()

    logger.info("Computing weighted majority vote.")
    pred_df = adata.obs[prediction_cols].astype(str)

    if prob_cols:
        weights_arr = (
            adata.obs[prob_cols].astype(np.float32) * weights_arr
        )  # (n_cells, n_methods)
    else:
        weights_arr = pd.DataFrame([weights_arr] * len(pred_df), index=pred_df.index)

    results = pred_df.apply(
        _weighted_majority, axis=1, weights=weights_arr, tie_label=par["tie_label"]
    )
    consensus_predictions, consensus_scores = zip(*results)

    logger.info("Writing output data.")
    adata.obs[par["output_obs_predictions"]] = pd.Categorical(consensus_predictions)
    adata.obs[par["output_obs_score"]] = consensus_scores

    logger.info("Writing output data...")
    write_h5ad_to_h5mu_with_compression(
        par["output"], par["input"], par["modality"], adata, par["output_compression"]
    )


if __name__ == "__main__":
    main()
