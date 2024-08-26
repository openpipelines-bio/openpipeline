import sys
import pytest
import mudata
import numpy as np

## VIASH START
meta = {
    'executable': './target/native/integrate/harmonypy/harmonypy',
    'resources_dir': './resources_test/pbmc_1k_protein_v3/'
}
## VIASH END

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3_mms.h5mu"

def test_harmonypy(run_component, tmp_path):
    output_path = tmp_path / "output.h5mu"

    # run component
    run_component([
        "--input", input_file,
        "--modality", "rna",
        "--obsm_input", "X_pca",
        "--obsm_output", "X_pca_int",
        "--obs_covariates", "harmony_integration_leiden_1.0",
        "--output", str(output_path),
        "--output_compression", "gzip"])
    assert output_path.is_file()

    # check output
    input_data = mudata.read_h5mu(input_file)
    output_data = mudata.read_h5mu(output_path)
    np.testing.assert_array_equal(output_data.mod['rna'].X.data, input_data.mod['rna'].X.data)
    np.testing.assert_array_equal(input_data.mod['rna'].obsm['X_pca'], output_data.mod['rna'].obsm['X_pca'])
    assert 'X_pca_int' in output_data.mod['rna'].obsm
    assert output_data.mod['rna'].obsm['X_pca_int'].shape == input_data.mod['rna'].obsm['X_pca'].shape

if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))