import muon as mu
import numpy as np

### VIASH START
par = {
  'input': 'resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu',
  'modality': ['rna'],
  'obs_filter': ['filter_none', 'filter_with_random'],
  'var_filter': ['filter_with_random'],
  'output': 'output.h5mu'
}

mdata = mu.read_h5mu(par["input"])
mdata.mod['rna'].obs["filter_none"] = np.repeat(True, mdata.mod['rna'].n_obs)
mdata.mod['rna'].obs["filter_with_random"] = np.random.choice(a=[False, True], size=mdata.mod['rna'].n_obs)
mdata.mod['rna'].var["filter_with_random"] = np.random.choice(a=[False, True], size=mdata.mod['rna'].n_vars)
mod = 'rna'
### VIASH END

print(f"Reading {par['input']}")
mdata = mu.read_h5mu(par["input"])

for mod in par["modality"]:
    print(f"Processing modality '{mod}'")

    obs_filt = np.repeat(True, mdata.mod[mod].n_obs)
    var_filt = np.repeat(True, mdata.mod[mod].n_vars)

    for obs_name in par["obs_filter"]:
        print(f"Filtering modality '{mod}' observations by .obs['{obs_name}']")
        if obs_name in mdata.mod[mod].obs:
            obs_filt &= mdata.mod[mod].obs[obs_name]
        else:
            print(f"Warning: .mod['{mod}'].obs['{obs_name}'] does not exist. Skipping.")

    for var_name in par["var_filter"]:
        print(f"Filtering modality '{mod}' variables by .var['{var_name}']")
        if var_name in mdata.mod[mod].var:
            var_filt &= mdata.mod[mod].var[var_name]
        else:
            print(f"Warning: .mod['{mod}'.var['{var_name}'] does not exist. Skipping.")
    
    mdata.mod[mod] = mdata.mod[mod][obs_filt, var_filt].copy()

print("Writing h5mu to file")
mdata.write_h5mu(par["output"])
