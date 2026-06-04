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
expected_mod = ["prot", "rna"]
assert all(key in list(input.mod.keys()) for key in expected_mod), (
    f"Query obs columns should be: {expected_mod}, found: {input.mod.keys()}."
)
assert input.shape[0] == output.shape[0], "Number of observations differ"
assert output.shape[1] < input.shape[1], "Features weren't filtered"

# Check rna modality
expected_rna_obs_keys = [
    "sample_id",
    "fraction_mitochondrial",
    "fraction_ribosomal",
    "filter_mitochondrial",
    "filter_ribosomal",
    "filter_with_counts",
    "filter_with_scrublet",
]
expected_rna_var_keys = ["mitochondrial", "ribosomal", "filter_with_counts"]
assert all(key in list(output.mod["rna"].obs) for key in expected_rna_obs_keys), (
    f"Query obs columns should be: {expected_rna_obs_keys}, found: {output.mod['rna'].obs.keys()}."
)
assert all(key in list(output.mod["rna"].var) for key in expected_rna_var_keys), (
    f"Query var columns should be: {expected_rna_var_keys}, found: {output.mod['rna'].var.keys()}."
)
assert input.mod["rna"].shape[0] == output.mod["rna"].shape[0], (
    "Number of observations in RNA modality differ"
)
assert output.mod["rna"].shape[1] < input.mod["rna"].shape[1], (
    "RNA features weren't filtered"
)

# Check prot modality
expected_prot_obs_keys = ["sample_id", "filter_with_counts"]
expected_prot_var_keys = ["filter_with_counts"]
assert all(key in list(output.mod["prot"].obs) for key in expected_prot_obs_keys), (
    f"Query obs columns should be: {expected_prot_obs_keys}, found: {output.mod['prot'].obs.keys()}."
)
assert all(key in list(output.mod["prot"].var) for key in expected_prot_var_keys), (
    f"Query var columns should be: {expected_prot_var_keys}, found: {output.mod['prot'].var.keys()}."
)
assert input.mod["prot"].shape == output.mod["prot"].shape, (
    "Protein modality shapes differ"
)

# Check that intersect_obs was applied: all modalities share the same observations
rna_obs = set(output.mod["rna"].obs_names)
prot_obs = set(output.mod["prot"].obs_names)
assert rna_obs == prot_obs, (
    f"After intersect_obs, RNA and prot modalities should share the same observations. "
    f"RNA-only: {len(rna_obs - prot_obs)}, prot-only: {len(prot_obs - rna_obs)}."
)

print("Test successful!", flush=True)
