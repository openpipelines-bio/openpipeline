from __future__ import annotations
import sys
import numpy as np
import pandas as pd
from mudata import read_h5mu

### VIASH START
par = {
    "input": "test_data/proportions_input.h5mu",
    "modality": "rna",
    "obs_participant_id": "participant_id",
    "obs_subpopulation": "subpopulation",
    "output": "test_data/proportions_output.h5mu",
    "uns_output": "proportions",
    "obsm_output": "proportions",
    "output_compression": None,
}
meta = {
    "resources_dir": ".",
}
### VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()


def main():
    logger.info("Reading input file '%s'.", par["input"])
    mdata = read_h5mu(par["input"])

    modality = par["modality"]
    if modality not in mdata.mod:
        raise ValueError(
            f"Modality '{modality}' not found in the MuData object. "
            f"Available modalities: {list(mdata.mod.keys())}"
        )
    adata = mdata.mod[modality]

    participant_col = par["obs_participant_id"]
    subpop_col = par["obs_subpopulation"]

    for col in (participant_col, subpop_col):
        if col not in adata.obs.columns:
            raise ValueError(
                f"Column '{col}' not found in .obs. "
                f"Available columns: {list(adata.obs.columns)}"
            )

    logger.info(
        "Computing proportions: '%s' × '%s'.",
        participant_col,
        subpop_col,
    )

    # Count cells per (participant, subpopulation)
    counts = (
        adata.obs.groupby([participant_col, subpop_col], observed=True)
        .size()
        .unstack(fill_value=0)
    )
    # Normalise rows to proportions (sum = 1 per participant)
    proportions = counts.div(counts.sum(axis=1), axis=0)

    n_participants = proportions.shape[0]
    n_subpops = proportions.shape[1]
    logger.info(
        "Proportion matrix shape: %d participants × %d subpopulations.",
        n_participants,
        n_subpops,
    )

    # Store in .uns as a serialisable dict
    uns_key = par["uns_output"]
    mdata.mod[modality].uns[uns_key] = proportions.to_dict()
    logger.info("Stored proportion matrix in .uns['%s'].", uns_key)

    # Store per-cell proportions in .obsm
    # Each cell gets the proportion row corresponding to its participant
    obsm_key = par["obsm_output"]
    participant_ids = adata.obs[participant_col].values
    obsm_matrix = np.array(
        [
            proportions.loc[pid].values
            if pid in proportions.index
            else np.zeros(n_subpops)
            for pid in participant_ids
        ],
        dtype=np.float64,
    )
    # Wrap in a DataFrame so column names (subpopulation labels) are preserved
    obsm_df = pd.DataFrame(
        obsm_matrix,
        index=adata.obs_names,
        columns=proportions.columns.astype(str),
    )
    mdata.mod[modality].obsm[obsm_key] = obsm_df
    logger.info(
        "Stored per-cell proportion vectors in .obsm['%s'] (shape: %s).",
        obsm_key,
        obsm_df.shape,
    )

    logger.info("Writing output to '%s'.", par["output"])
    mdata.write_h5mu(par["output"], compression=par["output_compression"])
    logger.info("Finished.")


if __name__ == "__main__":
    main()
