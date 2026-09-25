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


def _check_columns(obs, columns):
    """Fail when a requested `.obs` column is missing."""
    for col in columns:
        if col not in obs.columns:
            raise ValueError(
                f"Column '{col}' not found in .obs. "
                f"Available columns: {list(obs.columns)}"
            )


def _add_single_group(obs):
    """Add a constant group column, so all cells form one group named `all`."""
    group_col = "_obs_group"
    while group_col in obs.columns:
        group_col += "_"

    return obs.assign(**{group_col: "all"}), group_col


def _label_strata(obs, label_col, within_col):
    """Map every label to the one `--obs_normalize_within` value it occurs with."""
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

    return {str(label): str(values[0]) for label, values in label_to_within.items()}


def _normalise_within(counts, strata):
    """Divide each label's count by the group's total over the labels of its stratum."""
    denominators = pd.DataFrame(index=counts.index, columns=counts.columns)
    for stratum in sorted(set(strata.values())):
        members = [c for c in counts.columns if strata[str(c)] == stratum]
        totals = counts[members].sum(axis=1)
        for member in members:
            denominators[member] = totals

    return counts.div(denominators.astype(float)).fillna(0.0)


def compute_proportions(obs, group_col, label_col, within_col):
    """Group x label proportion matrix."""
    counts = (
        obs.groupby([group_col, label_col], observed=True).size().unstack(fill_value=0)
    )

    if within_col is None:
        proportions = counts.div(counts.sum(axis=1), axis=0)
    else:
        strata = _label_strata(obs, label_col, within_col)
        proportions = _normalise_within(counts, strata)
        logger.info(
            "Normalised within '%s' (%d strata); rows sum to the number of strata "
            "present in a group, not to 1.",
            within_col,
            len(set(strata.values())),
        )

    proportions.index = proportions.index.astype(str)
    proportions.columns = proportions.columns.astype(str)

    return proportions


def main():
    if not par["output"] and not par["output_csv"]:
        raise ValueError("Give at least one of --output or --output_csv.")

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
    _check_columns(adata.obs, required_cols)

    obs = adata.obs
    if group_col is None:
        obs, group_col = _add_single_group(obs)
        logger.info(
            "No --obs_group given; computing overall '%s' proportions over all cells.",
            label_col,
        )
    else:
        logger.info("Computing proportions: '%s' x '%s'.", group_col, label_col)

    proportions = compute_proportions(obs, group_col, label_col, within_col)
    logger.info("Proportion matrix shape: %d groups x %d labels.", *proportions.shape)

    if par["output_csv"]:
        csv_df = proportions.copy()
        csv_df.index.name = par["obs_group"] or "group"
        csv_df.to_csv(par["output_csv"])
        logger.info("Written proportion matrix to '%s'.", par["output_csv"])

    if par["output"]:
        # Store in .uns as a DataFrame (groups x labels), not a nested dict
        mdata.mod[modality].uns[par["uns_output"]] = proportions
        logger.info("Stored proportion matrix in .uns['%s'].", par["uns_output"])
        logger.info("Writing output to '%s'.", par["output"])
        mdata.write_h5mu(par["output"], compression=par["output_compression"])

    logger.info("Finished.")


if __name__ == "__main__":
    main()
