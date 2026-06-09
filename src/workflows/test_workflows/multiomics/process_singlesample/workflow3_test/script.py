import mudata as mu
import numpy as np

##VIASH START
par = {
    "input": "output.h5mu",
    "orig_input": "input.h5mu",
    "scrublet_score_threshold": 0.1,
}

meta = {"resources_dir": "resources_test/pbmc_1k_protein_v3"}

##VIASH END

print("Loading data", flush=True)
input = mu.read_h5mu(par["orig_input"])
output = mu.read_h5mu(par["input"])

assert input.n_mod == output.n_mod, "Number of modalities differ"
assert input.mod.keys() == output.mod.keys(), "Modalities differ"
expected_mod = ["prot", "rna"]
assert all(key in list(output.mod.keys()) for key in expected_mod), (
    f"Output modalities should be: {expected_mod}, found: {output.mod.keys()}."
)

rna_obs = output.mod["rna"].obs
assert "scrublet_doublet_score" in rna_obs.columns, (
    "scrublet_doublet_score column missing from rna .obs"
)
assert "filter_with_scrublet" in rna_obs.columns, (
    "filter_with_scrublet column missing from rna .obs"
)

# The manual threshold means: cells with doublet score > threshold are
# classified as doublets and are removed by the do_filter step that follows
# scrublet. Every cell remaining in the output must therefore be a non-doublet
# (scrublet_doublet_score <= threshold, filter_with_scrublet == True).
threshold = par["scrublet_score_threshold"]
doublet_scores = rna_obs["scrublet_doublet_score"].to_numpy(dtype=float)
keep_cells = rna_obs["filter_with_scrublet"].to_numpy(dtype=bool)

assert not np.isnan(doublet_scores).all(), (
    "scrublet_doublet_score column should not be all NaN when a manual threshold is set"
)

valid = ~np.isnan(doublet_scores)
assert (doublet_scores[valid] <= threshold).all(), (
    f"Cells with a doublet score above the manual threshold of {threshold} should "
    f"have been removed, but {(doublet_scores[valid] > threshold).sum()} remain in "
    f"the output."
)
assert keep_cells[valid].all(), (
    "All cells remaining after filtering should be tagged as kept "
    "(filter_with_scrublet == True)."
)

# Confirm that doublet removal actually reduced the number of RNA observations
# relative to the unfiltered input.
assert output.mod["rna"].n_obs < input.mod["rna"].n_obs, (
    f"Expected doublet removal to reduce the number of RNA observations below the "
    f"input ({input.mod['rna'].n_obs}), but got {output.mod['rna'].n_obs}."
)

print("Test successful!", flush=True)
