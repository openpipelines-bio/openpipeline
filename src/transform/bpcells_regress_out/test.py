import sys
import pytest
import mudata as mu
import numpy as np


orig_input = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"

def test_regress_out(run_component, tmp_path):
    input = tmp_path / "input.h5mu"
    output = tmp_path / "output.h5mu"

    # add column to obs
    mu_orig_input = mu.read_h5mu(orig_input)
    mu_orig_input.obs["var"] = np.random.rand(mu_orig_input.n_obs)
    mu_orig_input.mod["rna"].obs["var"] = mu_orig_input.obs["var"]
    mu_orig_input.mod["prot"].obs["var"] = mu_orig_input.obs["var"]
    mu_orig_input.write_h5mu(input)

    # execute command
    cmd_pars = [
        "--input", input,
        "--output", str(output),
        "--obs_keys", "var",
        "--output_compression", "gzip"
    ]
    run_component(cmd_pars)

    assert output.is_file(), "No output was created."

    mu_input = mu.read_h5mu(input)
    mu_output = mu.read_h5mu(output)

    assert "rna" in mu_output.mod, 'Output should contain data.mod["prot"].'
    assert "prot" in mu_output.mod, 'Output should contain data.mod["prot"].'

    rna_in = mu_input.mod["rna"]
    rna_out = mu_output.mod["rna"]
    prot_in = mu_input.mod["prot"]
    prot_out = mu_output.mod["prot"]

    assert rna_in.shape == rna_out.shape, "Should have same shape as before"
    assert prot_in.shape == prot_out.shape, "Should have same shape as before"

    assert np.mean(rna_in.X) != np.mean(rna_out.X), "Expression should have changed"

if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))