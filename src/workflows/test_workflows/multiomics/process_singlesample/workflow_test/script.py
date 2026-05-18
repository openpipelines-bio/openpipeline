import mudata as mu

##VIASH START
par = {"input": "output.h5mu", "orig_input": "input.5mu"}

meta = {"resources_dir": "resources_test/pbmc_1k_protein_v3"}

##VIASH END

print("Loading data", flush=True)
input = mu.read_h5mu(par["orig_input"])
output = mu.read_h5mu(par["input"])

assert input.n_mod == output.n_mod, "Number of modalities differ"
assert input.mod.keys() == output.mod.keys(), "Modalities differ"
expected_mod = ["atac", "rna"]
assert all(key in list(input.mod.keys()) for key in expected_mod), (
    f"Query obs columns should be: {expected_mod}, found: {input.mod.keys()}."
)
assert input.shape == output.shape, (
    f"Expected shape {input.shape} but got {output.shape}."
)

expected_rna_obs_keys = [
    "filter_with_counts",
    "scrublet_doublet_score",
    "filter_with_scrublet",
]

# Check rna modality
assert all(key in list(output.mod["rna"].obs) for key in expected_rna_obs_keys), (
    f"Query obs columns should be: {expected_rna_obs_keys}, found: {output.mod['rna'].obs.keys()}."
)
assert input.mod["rna"].shape == output.mod["rna"].shape, (
    f"Expected shape {input.mod['rna'].shape} but got {output.mod['rna'].shape}."
)
assert (input.mod["rna"].X != output.mod["rna"].X).nnz == 0, "RNA .X matrices differ"

# Check atac modality
assert not any(key in list(output.mod["atac"].obs) for key in expected_rna_obs_keys), (
    f"Query obs columns should NOT contain: {expected_rna_obs_keys}, but found these present: {[key for key in expected_rna_obs_keys if key in output.mod['atac'].obs.keys()]}."
)
assert input.mod["atac"].shape == output.mod["atac"].shape, (
    f"Expected shape {input.mod['atac'].shape} but got {output.mod['atac'].shape}."
)
assert (input.mod["atac"].X != output.mod["atac"].X).nnz == 0, "ATAC .X matrices differ"

print("Test successful!", flush=True)
