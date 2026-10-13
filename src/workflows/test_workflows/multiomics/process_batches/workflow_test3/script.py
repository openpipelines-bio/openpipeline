import mudata as mu

##VIASH START
par = {
    "input": "test.h5mu",
    "orig_input": "resources_test/10x_5k_anticmv/5k_human_antiCMV_T_TBNK_connect_mms.h5mu",
}

meta = {"resources_dir": "resources_test/10x_5k_anticmv"}

##VIASH END

print("Loading data", flush=True)
input = mu.read_h5mu(par["orig_input"])
output = mu.read_h5mu(par["input"])

assert input.n_mod == output.n_mod, "Number of modalities differ"
assert input.mod.keys() == output.mod.keys(), "Modalities differ"

# All modalities should share the same set of observations after intersect_obs
mod_obs_names = {name: set(mod.obs_names) for name, mod in output.mod.items()}
reference_name, reference_obs = next(iter(mod_obs_names.items()))
for name, obs in mod_obs_names.items():
    assert obs == reference_obs, (
        f"After intersect_obs, modality '{name}' should share the same observations "
        f"as '{reference_name}'. "
        f"{name}-only: {len(obs - reference_obs)}, "
        f"{reference_name}-only: {len(reference_obs - obs)}."
    )

# The intersected observations must be a subset of each original modality's observations
for name, mod in input.mod.items():
    assert reference_obs.issubset(set(mod.obs_names)), (
        f"Intersected observations should be a subset of the original '{name}' observations."
    )

print("Test successful!", flush=True)
