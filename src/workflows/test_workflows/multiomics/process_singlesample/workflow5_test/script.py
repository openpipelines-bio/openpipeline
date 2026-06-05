import mudata as mu

##VIASH START
par = {
    "input": "output.h5mu",
    "orig_input": "input.h5mu",
    "output_obs_num_nonzero_vars": "num_nonzero_vars",
    "output_obs_total_counts_vars": "total_counts",
    "output_var_num_nonzero_obs": "num_nonzero_obs",
    "output_var_total_counts_obs": "total_counts",
    "output_var_obs_mean": "obs_mean",
    "output_var_pct_dropout": "pct_dropout",
    "top_n_vars": [50, 100, 200, 500],
    "log1p_transform": True,
    "var_name_mitochondrial_genes": "mitochondrial",
    "var_name_ribosomal_genes": "ribosomal",
}

meta = {"resources_dir": "resources_test/pbmc_1k_protein_v3"}

##VIASH END

# Default slot names and top_n_vars, matching the process_singlesample defaults.
# These are the names that must NOT appear once a slot has been renamed.
DEFAULT_OBS_TOTAL_COUNTS = "total_counts"
DEFAULT_OBS_NUM_NONZERO = "num_nonzero_vars"
DEFAULT_VAR_OBS_MEAN = "obs_mean"
DEFAULT_VAR_TOTAL_COUNTS = "total_counts"
DEFAULT_VAR_NUM_NONZERO = "num_nonzero_obs"
DEFAULT_VAR_PCT_DROPOUT = "pct_dropout"
DEFAULT_TOP_N_VARS = [50, 100, 200, 500]

log1p = par["log1p_transform"]


def top_n_column(n):
    return f"pct_of_counts_in_top_{n}_vars"


print("Loading data", flush=True)
input = mu.read_h5mu(par["orig_input"])
output = mu.read_h5mu(par["input"])

assert input.n_mod == output.n_mod, "Number of modalities differ"
assert input.mod.keys() == output.mod.keys(), "Modalities differ"
expected_mod = ["prot", "rna"]
assert all(key in list(output.mod.keys()) for key in expected_mod), (
    f"Expected modalities {expected_mod}, found: {list(output.mod.keys())}."
)

# QC metric columns calculated for every processed modality.
common_obs_present = [
    par["output_obs_num_nonzero_vars"],
    par["output_obs_total_counts_vars"],
]
common_obs_present += [top_n_column(n) for n in par["top_n_vars"]]
if log1p:
    common_obs_present += [
        f"log1p_{par['output_obs_num_nonzero_vars']}",
        f"log1p_{par['output_obs_total_counts_vars']}",
    ]

common_var_present = [
    par["output_var_obs_mean"],
    par["output_var_total_counts_obs"],
    par["output_var_num_nonzero_obs"],
    par["output_var_pct_dropout"],
]
if log1p:
    common_var_present += [
        f"log1p_{par['output_var_obs_mean']}",
        f"log1p_{par['output_var_total_counts_obs']}",
    ]

# Per-cell proportions for boolean .var columns are only calculated for the rna
# modality (mitochondrial and ribosomal genes). These column names are fixed and
# not affected by the output slot arguments.
rna_obs_present = list(common_obs_present)
for metric in [par["var_name_mitochondrial_genes"], par["var_name_ribosomal_genes"]]:
    rna_obs_present += [f"total_counts_{metric}", f"pct_{metric}"]
    if log1p:
        rna_obs_present.append(f"log1p_total_counts_{metric}")

# Default-named columns that must be absent when the corresponding slot is renamed.
absent_obs = []
if par["output_obs_num_nonzero_vars"] != DEFAULT_OBS_NUM_NONZERO:
    absent_obs.append(DEFAULT_OBS_NUM_NONZERO)
    if log1p:
        absent_obs.append(f"log1p_{DEFAULT_OBS_NUM_NONZERO}")
if par["output_obs_total_counts_vars"] != DEFAULT_OBS_TOTAL_COUNTS:
    absent_obs.append(DEFAULT_OBS_TOTAL_COUNTS)
    if log1p:
        absent_obs.append(f"log1p_{DEFAULT_OBS_TOTAL_COUNTS}")
absent_obs += [
    top_n_column(n) for n in DEFAULT_TOP_N_VARS if n not in par["top_n_vars"]
]

absent_var = []
if par["output_var_obs_mean"] != DEFAULT_VAR_OBS_MEAN:
    absent_var.append(DEFAULT_VAR_OBS_MEAN)
    if log1p:
        absent_var.append(f"log1p_{DEFAULT_VAR_OBS_MEAN}")
if par["output_var_total_counts_obs"] != DEFAULT_VAR_TOTAL_COUNTS:
    absent_var.append(DEFAULT_VAR_TOTAL_COUNTS)
    if log1p:
        absent_var.append(f"log1p_{DEFAULT_VAR_TOTAL_COUNTS}")
if par["output_var_num_nonzero_obs"] != DEFAULT_VAR_NUM_NONZERO:
    absent_var.append(DEFAULT_VAR_NUM_NONZERO)
if par["output_var_pct_dropout"] != DEFAULT_VAR_PCT_DROPOUT:
    absent_var.append(DEFAULT_VAR_PCT_DROPOUT)

expected = {
    "rna": {"obs": rna_obs_present, "var": common_var_present},
    "prot": {"obs": common_obs_present, "var": common_var_present},
}

for modality, present in expected.items():
    mod_obs = list(output.mod[modality].obs)
    mod_var = list(output.mod[modality].var)

    missing_obs = [k for k in present["obs"] if k not in mod_obs]
    missing_var = [k for k in present["var"] if k not in mod_var]
    assert not missing_obs, (
        f"Missing QC metric obs columns for {modality}: {missing_obs}. Found: {mod_obs}."
    )
    assert not missing_var, (
        f"Missing QC metric var columns for {modality}: {missing_var}. Found: {mod_var}."
    )

    # When non-default slots are used, the default-named columns must be absent.
    unexpected_obs = [k for k in absent_obs if k in mod_obs]
    unexpected_var = [k for k in absent_var if k in mod_var]
    assert not unexpected_obs, (
        f"Default-named obs columns should be absent for {modality} when slots are "
        f"renamed, but found: {unexpected_obs}."
    )
    assert not unexpected_var, (
        f"Default-named var columns should be absent for {modality} when slots are "
        f"renamed, but found: {unexpected_var}."
    )

print("Test successful!", flush=True)
