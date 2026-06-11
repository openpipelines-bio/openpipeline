import sys
import mudata as mu
import pandas as pd

### VIASH START
par = {
    "input": "input.h5mu",
    "modality": "rna",
    "prot_modality": "prot",
    "sample_id_column": "sample_id",
    "rna_filter_columns": [
        "filter_counts_rna",
        "filter_quantile_rna",
        "filter_mito_rna",
    ],
    "scrublet_filter_column": "filter_scrublet",
    "prot_filter_columns": ["filter_counts_prot", "filter_quantile_prot"],
    "output": "cell_counts.tsv",
    "output_column_rna": "cell_count_after_rna_filter",
    "output_column_scrublet": "cell_count_after_scrublet_filter",
    "output_column_prot": "cell_count_after_protein_filter",
    "output_column_all": "cell_count_after_all_filter",
}
### VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()


def combine_keep_columns(obs, columns):
    """Logical AND of boolean keep-columns. Missing values are treated as not-kept."""
    keep = pd.Series(True, index=obs.index)
    for column in columns:
        if column not in obs:
            raise ValueError(f"Column '{column}' was not found in .obs.")
        keep &= obs[column].fillna(False).astype(bool)
    return keep


logger.info("Reading modality '%s' from %s", par["modality"], par["input"])
rna_obs = mu.read_h5ad(par["input"], mod=par["modality"]).obs

# Determine the per-sample grouping. When the sample id column is absent, treat the whole
# input as a single sample.
if par["sample_id_column"] in rna_obs:
    groups = rna_obs[par["sample_id_column"]]
else:
    logger.info(
        "Column '%s' not found in .obs, treating all observations as a single sample.",
        par["sample_id_column"],
    )
    groups = pd.Series("sample", index=rna_obs.index)

# RNA keep-mask (AND of the RNA filter columns).
rna_keep = combine_keep_columns(rna_obs, par["rna_filter_columns"] or [])

# Scrublet keep-mask (optional).
report_scrublet = par["scrublet_filter_column"] is not None
if report_scrublet:
    scrublet_keep = combine_keep_columns(rna_obs, [par["scrublet_filter_column"]])
else:
    scrublet_keep = pd.Series(True, index=rna_obs.index)

# Protein keep-mask (optional), aligned to the RNA observations.
report_prot = par["prot_modality"] is not None
if report_prot:
    logger.info("Reading modality '%s' from %s", par["prot_modality"], par["input"])
    prot_obs = mu.read_h5ad(par["input"], mod=par["prot_modality"]).obs
    prot_keep = combine_keep_columns(prot_obs, par["prot_filter_columns"] or [])
    # Observations not present in the protein modality cannot pass the protein filter.
    prot_keep = prot_keep.reindex(rna_obs.index, fill_value=False)
else:
    prot_keep = pd.Series(True, index=rna_obs.index)

# The combined mask: cells passing every filter (the intersection).
all_keep = rna_keep & scrublet_keep & prot_keep

logger.info("Aggregating cell counts per sample")
counts = pd.DataFrame({par["sample_id_column"]: groups})
counts[par["output_column_rna"]] = rna_keep
if report_scrublet:
    counts[par["output_column_scrublet"]] = scrublet_keep
if report_prot:
    counts[par["output_column_prot"]] = prot_keep
counts[par["output_column_all"]] = all_keep

# Summing booleans per group yields the number of cells passing each stage.
report = (
    counts.groupby(par["sample_id_column"], sort=True).sum().astype(int).reset_index()
)

logger.info("Writing cell count report to %s", par["output"])
report.to_csv(par["output"], sep="\t", index=False)
