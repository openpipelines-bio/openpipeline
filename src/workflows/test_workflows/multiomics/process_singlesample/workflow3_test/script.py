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
# classified as doublets (filter_with_scrublet == False), the rest are kept
# (filter_with_scrublet == True).
threshold = par["scrublet_score_threshold"]
doublet_scores = rna_obs["scrublet_doublet_score"].to_numpy(dtype=float)
keep_cells = rna_obs["filter_with_scrublet"].to_numpy(dtype=bool)

assert not np.isnan(doublet_scores).all(), (
    "scrublet_doublet_score column should not be all NaN when a manual threshold is set"
)

valid = ~np.isnan(doublet_scores)
expected_keep = doublet_scores[valid] <= threshold
actual_keep = keep_cells[valid]
assert np.array_equal(actual_keep, expected_keep), (
    f"filter_with_scrublet should match (scrublet_doublet_score <= {threshold}). "
    f"Mismatches: {(actual_keep != expected_keep).sum()} out of {valid.sum()} cells."
)

assert (doublet_scores[valid] > threshold).any(), (
    f"Expected at least one cell with doublet score above the manual threshold "
    f"of {threshold} to confirm the threshold was applied."
)

print("Test successful!", flush=True)
