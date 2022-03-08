from os import path
import subprocess
import muon

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


# TEST 1: filtering a little bit
print("> Check filtering a little bit")
out = subprocess.check_output([
        f"./{meta['functionality_name']}", 
        "--input", input_path, 
        "--output", "output-1.h5mu",
        "--min_cells_per_gene", "3"
]).decode("utf-8")

assert path.exists("output-1.h5mu"), "Output file not found"
mu_out = muon.read_h5mu("output-1.h5mu")

assert "filter_with_counts" in mu_out.mod["rna"].obs
assert "filter_with_counts" in mu_out.mod["rna"].var
new_obs = mu_out.mod['rna'].n_obs
new_vars = mu_out.mod['rna'].n_vars
assert new_obs == orig_obs, "No RNA obs should have been filtered"
assert new_vars == orig_vars, "No RNA vars should have been filtered"
assert mu_out.mod['prot'].n_obs == orig_prot_obs, "No prot obs should have been filtered"
assert mu_out.mod['prot'].n_vars == orig_prot_vars, "No prot vars should have been filtered"

assert list(mu_in.mod['rna'].var['feature_types'].cat.categories) == ["Gene Expression"], "Feature types of RNA modality should be Gene Expression"
assert list(mu_in.mod['prot'].var['feature_types'].cat.categories) == ["Antibody Capture"], "Feature types of prot modality should be Antibody Capture"


# TEST 2: filering a lot
print("> Check filtering a lot")
out = subprocess.check_output([
        f"./{meta['functionality_name']}", 
        "--input", input_path, 
        "--output", "output-2.h5mu",
        "--modality", "rna:prot",
        "--min_cells_per_gene", "100",
        "--min_counts", "200", 
        "--max_counts", "5000000",
        "--min_genes_per_cell", "200", 
        "--max_genes_per_cell", "1500000", 
        "--min_cells_per_gene", "10",
        "--min_fraction_mito", "0",
        "--max_fraction_mito", "0.2",
        "--do_subset"
]).decode("utf-8")

assert path.exists("output-2.h5mu"), "Output file not found"
mu_out = muon.read_h5mu("output-2.h5mu")

new_new_obs = mu_out.mod['rna'].n_obs
new_new_vars = mu_out.mod['rna'].n_vars
assert new_new_obs <= new_obs, "More cells should have been filtered"
assert new_new_vars <= new_vars, "More genes should have been filtered"
assert mu_out.mod['prot'].n_obs <= orig_prot_obs, "Some prot obs should have been filtered"
assert mu_out.mod['prot'].n_vars <= orig_prot_vars, "Some prot vars should have been filtered"

assert list(mu_in.mod['rna'].var['feature_types'].cat.categories) == ["Gene Expression"], "Feature types of RNA modality should be Gene Expression"
assert list(mu_in.mod['prot'].var['feature_types'].cat.categories) == ["Antibody Capture"], "Feature types of prot modality should be Antibody Capture"

# test succeeded
print("> All tests succeeded!")