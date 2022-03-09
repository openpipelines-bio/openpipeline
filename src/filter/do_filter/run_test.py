from os import path
import subprocess
import muon
import numpy as np

## VIASH START
meta = {
    'functionality_name': 'foo',
    'resources_dir': 'resources_test/'
}
## VIASH END

print("> Reading input file")
input_path = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"

mu_in = muon.read_h5mu(input_path)
orig_obs = mu_in.mod['rna'].n_obs
orig_vars = mu_in.mod['rna'].n_vars
orig_prot_obs = mu_in.mod['prot'].n_obs
orig_prot_vars = mu_in.mod['prot'].n_vars

ad_rna = mu_in.mod['rna']
print(f"  input: {ad_rna}")
ad_rna.obs["filter_none"] = np.repeat(True, ad_rna.n_obs)
ad_rna.obs["filter_with_random"] = np.random.choice(a=[False, True], size=ad_rna.n_obs)
ad_rna.var["filter_with_random"] = np.random.choice(a=[False, True], size=ad_rna.n_vars)

mu_in.write_h5mu("input_with_extra_columns.h5mu")

# TEST 1: filtering a little bit
print("> Check filtering a little bit")
out = subprocess.check_output([
        f"./{meta['functionality_name']}", 
        "--input", "input_with_extra_columns.h5mu", 
        "--output", "output-1.h5mu",
        "--obs_filter", "filter_none:filter_with_random",
        "--var_filter", "filter_with_random"
]).decode("utf-8")

assert path.exists("output-1.h5mu"), "Output file not found"
mu_out = muon.read_h5mu("output-1.h5mu")

print(f"  output1: {mu_out.mod['rna']}")
new_obs = mu_out.mod['rna'].n_obs
new_vars = mu_out.mod['rna'].n_vars
assert new_obs < orig_obs, "Some RNA obs should have been filtered"
assert new_vars < orig_vars, "Some RNA vars should have been filtered"


# TEST 2: filtering nothing
print("> Check filtering a little bit")
out = subprocess.check_output([
        f"./{meta['functionality_name']}", 
        "--input", "input_with_extra_columns.h5mu", 
        "--output", "output-2.h5mu",
        "--obs_filter", "filter_none"
]).decode("utf-8")

assert path.exists("output-2.h5mu"), "Output file not found"
mu_out = muon.read_h5mu("output-2.h5mu")

print(f"  output2: {mu_out.mod['rna']}")
new_obs = mu_out.mod['rna'].n_obs
new_vars = mu_out.mod['rna'].n_vars
assert new_obs == orig_obs, "No RNA obs should have been filtered"
assert new_vars == orig_vars, "No RNA vars should have been filtered"

# test succeeded
print("> All tests succeeded!")