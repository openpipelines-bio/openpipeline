from os import path
import subprocess
import muon
import logging
from sys import stdout

## VIASH START
meta = {
    'functionality_name': 'foo',
    'resources_dir': 'resources_test/'
}
## VIASH END

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

logger.info("> Reading input file")
input_path = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"

mu_in = muon.read_h5mu(input_path)
orig_obs = mu_in.mod['rna'].n_obs
orig_vars = mu_in.mod['rna'].n_vars
orig_prot_obs = mu_in.mod['prot'].n_obs
orig_prot_vars = mu_in.mod['prot'].n_vars


# TEST 1: filtering a little bit
logger.info("> Check filtering a little bit")
out = subprocess.check_output([
        f"./{meta['functionality_name']}", 
        "--input", input_path, 
        "--output", "output-1.h5mu",
        "--min_counts", "3"
]).decode("utf-8")

assert path.exists("output-1.h5mu"), "Output file not found"
mu_out = muon.read_h5mu("output-1.h5mu")

assert "filter_with_scrublet" in mu_out.mod["rna"].obs
new_obs = mu_out.mod['rna'].n_obs
new_vars = mu_out.mod['rna'].n_vars
assert new_obs == orig_obs, "No RNA obs should have been filtered"
assert new_vars == orig_vars, "No RNA vars should have been filtered"
assert mu_out.mod['prot'].n_obs == orig_prot_obs, "No prot obs should have been filtered"
assert mu_out.mod['prot'].n_vars == orig_prot_vars, "No prot vars should have been filtered"

assert list(mu_in.mod['rna'].var['feature_types'].cat.categories) == ["Gene Expression"], "Feature types of RNA modality should be Gene Expression"
assert list(mu_in.mod['prot'].var['feature_types'].cat.categories) == ["Antibody Capture"], "Feature types of prot modality should be Antibody Capture"


# TEST 2: filering a lot
logger.info("> Check filtering a lot")
out = subprocess.check_output([
        f"./{meta['functionality_name']}", 
        "--input", input_path, 
        "--output", "output-2.h5mu",
        "--modality", "rna",
        "--min_counts", "10",
        "--num_pca_components", "10",
        "--do_subset"
]).decode("utf-8")

assert path.exists("output-2.h5mu"), "Output file not found"
mu_out = muon.read_h5mu("output-2.h5mu")

new_new_obs = mu_out.mod['rna'].n_obs
new_new_vars = mu_out.mod['rna'].n_vars
assert new_new_obs <= new_obs, "More cells should have been filtered"
assert new_new_vars <= new_vars, "More genes should have been filtered"
assert mu_out.mod['prot'].n_obs == orig_prot_obs, "No prot obs should have been filtered"
assert mu_out.mod['prot'].n_vars == orig_prot_vars, "No prot vars should have been filtered"

assert list(mu_in.mod['rna'].var['feature_types'].cat.categories) == ["Gene Expression"], "Feature types of RNA modality should be Gene Expression"
assert list(mu_in.mod['prot'].var['feature_types'].cat.categories) == ["Antibody Capture"], "Feature types of prot modality should be Antibody Capture"

# test succeeded
logger.info("> All tests succeeded!")