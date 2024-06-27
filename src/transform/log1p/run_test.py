from os import path
import mudata as mu
import numpy as np
import scanpy as sc
import pandas as pd
import sys
import pytest
import sys
import uuid
from operator import attrgetter

## VIASH START
meta = {
    'functionality_name': 'lognorm',
    'resources_dir': 'resources_test/',
    'config': './src/transform/log1p/config.vsh.yaml',
    'executable': "../../executable/docker/transform/log1p/log1p"
}


## VIASH END

@pytest.fixture
def input_data():
    return mu.read_h5mu(f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu").copy()

@pytest.fixture
def random_h5mu_path(tmp_path):
    def wrapper():
        unique_filename = f"{str(uuid.uuid4())}.h5mu"
        temp_file = tmp_path / unique_filename
        return temp_file
    return wrapper

@pytest.mark.parametrize("output_layer", [None, "log_normalized"])
@pytest.mark.parametrize("input_layer", [None, "normalized"])
def test_1logp(run_component, input_data, output_layer, input_layer, random_h5mu_path):
    output = random_h5mu_path()
    if input_layer:
        mod = input_data.mod["rna"]
        mod.layers[input_layer] = mod.X.copy()
        # Overwrite the original layer to make sure
        # it is not accidentally used as input layer.
        mod.X[:] = 0
    input_path = random_h5mu_path()
    input_data.write(input_path)
    run_args = [
        "--input", input_path,
        "--output", output,
        "--output_compresion", "gzip"
    ]
    if output_layer:
        run_args.extend(["--output_layer", output_layer])
    if input_layer:
        run_args.extend(["--input_layer", input_layer])
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
    input_layer_data = rna_in.X if not input_layer else rna_in.layers[input_layer]
    assert np.mean(input_layer_data) != np.mean(get_output_layer(rna_out)), "Expression should have changed"

    nz_row, nz_col = input_layer_data.nonzero()
    row_corr = np.corrcoef(input_layer_data[nz_row[0],:].toarray().flatten(), 
                           get_output_layer(rna_out)[nz_row[0],:].toarray().flatten())[0,1]
    col_corr = np.corrcoef(input_layer_data[:,nz_col[0]].toarray().flatten(), 
                           get_output_layer(rna_out)[:,nz_col[0]].toarray().flatten())[0,1]
    assert row_corr > .1
    assert col_corr > .1

    assert 'log1p' in rna_out.uns

    # Make sure that the original input layer has not been overwritten
    layers_to_test = [None] + list(rna_in.layers.keys())
    for layer in layers_to_test:
        if layer != output_layer:
            in_data = sc.get.var_df(rna_in,
                                    keys=rna_in.obs_names.to_list(),
                                    layer=layer)
            out_data =  sc.get.var_df(rna_out,
                                      keys=rna_in.obs_names.to_list(),
                                      layer=layer)
            pd.testing.assert_frame_equal(in_data, out_data)
    

if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))