import sys
import subprocess
import pytest
from os import path
import mudata as mu
import numpy as np

## VIASH START
meta = {
    'name': 'lognorm',
    'resources_dir': 'resources_test/'
}
## VIASH END

input = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"

def test_run(run_component, tmp_path):
    output = tmp_path / "output.h5mu"
    cmd_pars = [
        "--input", input,
        "--output", str(output),
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

    nz_row, nz_col = rna_in.X.nonzero()
    row_corr = np.corrcoef(rna_in.X[nz_row[0],:].toarray().flatten(), rna_out.X[nz_row[0],:].toarray().flatten())[0,1]
    col_corr = np.corrcoef(rna_in.X[:,nz_col[0]].toarray().flatten(), rna_out.X[:,nz_col[0]].toarray().flatten())[0,1]
    assert row_corr > .1
    assert col_corr > .1
    
def test_target_sum(run_component, tmp_path):
    output = tmp_path / "output.h5mu"
    cmd_pars = [
        "--input", input,
        "--output", str(output),
        "--output_compression", "gzip",
        "--target_sum", "10000"
    ]
    run_component(cmd_pars)

    assert output.is_file(), "No output was created."

    mu_output = mu.read_h5mu(output)
    rna_out = mu_output.mod["rna"]

    assert np.all(np.abs(rna_out.X.sum(axis=1) - 10000) < 1), "Expression should have changed"

if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))