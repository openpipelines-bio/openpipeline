from __future__ import annotations
import sys
import pandas as pd
from mudata import read_h5mu

### VIASH START
par = {
    "input": "proportions_input.h5mu",
    "modality": "rna",
    "obs_group": "participant_id",
    "obs_label": "subpopulation",  # generic: any .obs label column
    "obs_normalize_within": None,
    "output": "proportions_output.h5mu",
    "output_csv": None,
    "uns_output": "proportions",
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

    within_col = par["obs_normalize_within"]

    required_cols = [label_col] if group_col is None else [group_col, label_col]
    if within_col is not None:
        required_cols.append(within_col)
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

    if within_col is None:
        # Normalise rows to proportions (sum = 1 per group)
        proportions = counts.div(counts.sum(axis=1), axis=0)
    else:
        # Normalise within each (group, --obs_normalize_within) stratum. Every label
        # belongs to exactly one stratum, so the denominator of a label's column is the
        # group's total over the labels sharing its stratum.
        label_to_within = (
            obs[[label_col, within_col]]
            .drop_duplicates()
            .groupby(label_col, observed=True)[within_col]
            .agg(lambda values: sorted(set(values)))
        )
        ambiguous = {
            str(label): values
            for label, values in label_to_within.items()
            if len(values) > 1
        }
        if ambiguous:
            raise ValueError(
                f"Every '{label_col}' value must belong to exactly one "
                f"'{within_col}' value, but these span several: {ambiguous}. "
                f"'--obs_normalize_within' expects a strict hierarchy."
            )
        strata = {
            str(label): str(values[0]) for label, values in label_to_within.items()
        }
        denominators = pd.DataFrame(index=counts.index, columns=counts.columns)
        for stratum in sorted(set(strata.values())):
            members = [c for c in counts.columns if strata[str(c)] == stratum]
            totals = counts[members].sum(axis=1)
            for member in members:
                denominators[member] = totals
        denominators = denominators.astype(float)
        proportions = counts.div(denominators).fillna(0.0)
        logger.info(
            "Normalised within '%s' (%d strata); rows sum to the number of strata "
            "present in a group, not to 1.",
            within_col,
            len(set(strata.values())),
        )

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
