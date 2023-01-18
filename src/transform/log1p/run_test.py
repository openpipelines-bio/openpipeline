import subprocess
from os import path
import mudata as mu
import numpy as np
import logging
import pytest
import sys
from operator import getitem, attrgetter


## VIASH START
meta = {
    'functionality_name': 'lognorm',
    'resources_dir': 'resources_test/',
    'config': '/home/di/code/openpipeline/src/transform/log1p/config.vsh.yaml',
    'executable': "../../target/docker/transform/log1p/log1p"
}


## VIASH END

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(sys.stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

@pytest.fixture
def input_path():
    return f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"

@pytest.mark.parametrize("output_layer", [None, "log_normalized"])
def test_1logp(run_component, input_path, output_layer):
    output = "output.h5mu"
    run_args = [
        f"./{meta['functionality_name']}",
        "--input", input_path,
        "--output", output
        ]
    if output_layer:
        run_args.extend(["--output_layer", output_layer])
    run_component(run_args)
    get_output_layer = attrgetter("X") if not output_layer else lambda x: getattr(x, 'layers')[output_layer]

    assert path.exists(output), "No output was created."

    mu_input = mu.read_h5mu(input_path)
    mu_output = mu.read_h5mu(output)

    assert "rna" in mu_output.mod, 'Output should contain data.mod["prot"].'
    assert "prot" in mu_output.mod, 'Output should contain data.mod["prot"].'

    rna_in = mu_input.mod["rna"]
    rna_out = mu_output.mod["rna"]
    prot_in = mu_input.mod["prot"]
    prot_out = mu_output.mod["prot"]

    assert rna_in.shape == rna_out.shape, "Should have same shape as before"
    assert prot_in.shape == prot_out.shape, "Should have same shape as before"

    assert np.mean(rna_in.X) != np.mean(get_output_layer(rna_out)), "Expression should have changed"

    nz_row, nz_col = rna_in.X.nonzero()
    row_corr = np.corrcoef(rna_in.X[nz_row[0],:].toarray().flatten(), get_output_layer(rna_out)[nz_row[0],:].toarray().flatten())[0,1]
    col_corr = np.corrcoef(rna_in.X[:,nz_col[0]].toarray().flatten(), get_output_layer(rna_out)[:,nz_col[0]].toarray().flatten())[0,1]
    assert row_corr > .1
    assert col_corr > .1

    assert 'log1p' in rna_out.uns

if __name__ == '__main__':
    sys.exit(pytest.main([__file__], plugins=["viashpy"]))