import sys
import pytest
from mudata import read_h5mu

## VIASH START
meta = {
    'executable': './target/docker/graph/bbknn/bbknn',
    'resources_dir': './resources_test/pbmc_1k_protein_v3/'
}
## VIASH END

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3_mms.h5mu"

def test_simple_integration(run_component, tmp_path):
    output_path = tmp_path / "output.h5mu"
    tmp_input_path = tmp_path / "input.h5mu"

    # create input data
    input_data = read_h5mu(input_file)
    mod = input_data.mod['rna']
    if 'connectivities' in mod.obsp:
        del mod.obsp['connectivities']
    if 'distances' in mod.obsp:
        del mod.obsp['distances']
    if 'neighbors' in mod.uns:
        del mod.uns['neighbors']
    input_data.write(tmp_input_path)

    # run component
    run_component([
        "--input", str(tmp_input_path),
        "--output", str(output_path),
        "--obs_batch", "leiden",
        "--obsm_input", "X_pca",
        "--output_compression", "gzip"
    ])
    assert output_path.exists()
    data = read_h5mu(output_path).mod['rna']
    assert "connectivities" in data.obsp
    assert "distances" in data.obsp
    assert "neighbors" in data.uns

if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))