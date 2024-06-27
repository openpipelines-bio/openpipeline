import sys
import pytest
import subprocess
import mudata as mu

## VIASH START
meta = {
    'functionality_name': 'foo',
    'resources_dir': '/resources_test/pbmc_1k_protein_v3/'
}
## VIASH END

input = meta["resources_dir"] + "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"

def test_umap(run_component, tmp_path):
    output = tmp_path / "output.h5mu"
    run_component([
        "--input", input,
        "--output",  str(output),
        "--obsm_output", "X_foo",
        "--num_components", "26",
        "--output_compression", "gzip"
    ])
    
    assert output.is_file(), "No output was created."
    data = mu.read_h5mu(output)

    # check whether umap was found
    assert "X_foo" in data.mod["rna"].obsm, "Check whether output was found in .obsm"
    assert data.mod["rna"].obsm["X_foo"].shape == (data.n_obs, 26), "Check whether output has correct shape"

def test_raise_if_uns_neighbor_is_missing(run_component, tmp_path):
    output = tmp_path / "output.h5mu"
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component([
            "--input", input,
            "--output", str(output),
            "--obsm_output", "X_foo",
            "--num_components", "26",
            "--uns_neighbors", "does_not_exist"
        ])
    assert not output.is_file(), "No output should be created."
    assert "ValueError: 'does_not_exist' was not found in .mod['rna'].uns." in \
        err.value.stdout.decode('utf-8')

if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))