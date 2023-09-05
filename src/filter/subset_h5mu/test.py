import sys
import pytest
import mudata as mu

## VIASH START
meta = {
    'functionality_name': './target/docker/filter/subset_h5mu/subset_h5mu',
    'resources_dir': 'resources_test/pbmc_1k_protein_v3/'
}
## VIASH END

input_path = f"{meta['resources_dir']}/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"

def test_filter_nothing(run_component, tmp_path):
    output_path = tmp_path / "output.h5mu"
    
    # run component
    run_component([
        "--input", input_path,
        "--output", str(output_path),
        "--number_of_observations", "100",
        "--output_compression", "gzip"
    ])
    
    assert output_path.is_file(), "Output file not found"

    # check output file
    mu_in = mu.read_h5mu(input_path)
    mu_out = mu.read_h5mu(output_path)

    orig_vars = mu_in.mod['rna'].n_vars
    orig_prot_obs = mu_in.mod['prot'].n_obs
    orig_prot_vars = mu_in.mod['prot'].n_vars

    new_obs = mu_out.mod['rna'].n_obs
    new_vars = mu_out.mod['rna'].n_vars
    
    assert new_obs == 100, "Output should only contain 100 observations"
    assert new_vars == orig_vars, "No RNA vars should have been filtered"
    assert mu_out.mod['prot'].n_obs == orig_prot_obs, "No prot obs should have been filtered"
    assert mu_out.mod['prot'].n_vars == orig_prot_vars, "No prot vars should have been filtered"
    assert list(mu_out.mod['rna'].var['feature_types'].cat.categories) == ["Gene Expression"], \
        "Feature types of RNA modality should be Gene Expression"
    assert list(mu_out.mod['prot'].var['feature_types'].cat.categories) == ["Antibody Capture"], \
        "Feature types of prot modality should be Antibody Capture"

if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
