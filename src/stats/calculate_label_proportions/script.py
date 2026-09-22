from __future__ import annotations
import sys
import numpy as np
import pandas as pd
from mudata import read_h5mu

### VIASH START
par = {
    "input": "proportions_input.h5mu",
    "modality": "rna",
    "obs_group": "participant_id",
    "obs_label": "subpopulation",  # generic: any .obs label column
    "output": "proportions_output.h5mu",
    "output_csv": None,
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

    group_col = par["obs_group"]
    label_col = par["obs_label"]

    required_cols = [label_col] if group_col is None else [group_col, label_col]
    for col in required_cols:
        if col not in adata.obs.columns:
            raise ValueError(
                f"Column '{col}' not found in .obs. "
                f"Available columns: {list(adata.obs.columns)}"
            )

    obs = adata.obs
    if group_col is None:
        # No grouping column: all cells form a single group, i.e. overall proportions.
        group_col = "_obs_group"
        while group_col in obs.columns:
            group_col += "_"
        obs = obs.assign(**{group_col: "all"})
        logger.info(
            "No --obs_group given; computing overall '%s' proportions over all cells.",
            label_col,
        )
    else:
        logger.info(
            "Computing proportions: '%s' x '%s'.",
            group_col,
            label_col,
        )

    # Count cells per (group, label)
    counts = (
        obs.groupby([group_col, label_col], observed=True).size().unstack(fill_value=0)
    )
    # Normalise rows to proportions (sum = 1 per group)
    proportions = counts.div(counts.sum(axis=1), axis=0)

    n_groups = proportions.shape[0]
    n_labels = proportions.shape[1]
    logger.info(
        "Proportion matrix shape: %d groups x %d labels.",
        n_groups,
        n_labels,
    )

    # Store in .uns as a DataFrame (groups x labels), not a nested dict
    uns_key = par["uns_output"]
    proportions.index = proportions.index.astype(str)
    proportions.columns = proportions.columns.astype(str)
    mdata.mod[modality].uns[uns_key] = proportions
    logger.info("Stored proportion matrix in .uns['%s'].", uns_key)

    # Optional per-cell copy in .obsm. Every cell of a group carries the same row, so this
    # is redundant by construction; it exists only for cell-level components that take an
    # .obsm matrix (dimred/phate). Not written unless --obsm_output is given.
    obsm_key = par["obsm_output"]
    if obsm_key:
        group_ids = obs[group_col].values
        obsm_matrix = np.array(
            [
                proportions.loc[str(gid)].values
                if str(gid) in proportions.index
                else np.zeros(n_labels)
                for gid in group_ids
            ],
            dtype=np.float64,
        )
        # Wrap in a DataFrame so column names (labels) are preserved
        obsm_df = pd.DataFrame(
            obsm_matrix,
            index=adata.obs_names,
            columns=proportions.columns,
        )
        mdata.mod[modality].obsm[obsm_key] = obsm_df
        logger.info(
            "Stored per-cell proportion vectors in .obsm['%s'] (shape: %s).",
            obsm_key,
            obsm_df.shape,
        )

    # Optional tabular copy for downstream components that do not read MuData
    if par["output_csv"]:
        csv_df = proportions.copy()
        csv_df.index.name = par["obs_group"] or "group"
        csv_df.to_csv(par["output_csv"])
        logger.info("Written proportion matrix to '%s'.", par["output_csv"])

    logger.info("Writing output to '%s'.", par["output"])
    mdata.write_h5mu(par["output"], compression=par["output_compression"])
    logger.info("Finished.")


if __name__ == "__main__":
    main()
