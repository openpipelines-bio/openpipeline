import sys
import pytest
import mudata as mu

## VIASH START
meta = {
    'resources_dir': 'resources_test/pbmc_1k_protein_v3/'
}
## VIASH END

input_path = f"{meta['resources_dir']}/pbmc_1k_protein_v3_mms.h5mu"

def test_subset_obsp(run_component, tmp_path):
    output_path = tmp_path / "output.h5mu"

    # run component
    run_component([
        "--input", input_path,
        "--output", str(output_path),
        "--input_obsp_key", "distances",
        "--input_obs_key", "leiden",
        "--input_obs_value", "1",
        "--output_obsm_key", "leiden_1"
    ])

    assert output_path.is_file(), "Output file not found"

    # check output file
    mu_out = mu.read_h5mu(output_path)

    assert "leiden_1" in mu_out.mod["rna"].obsm, "Output should contain leiden_1 in .obsm"
    assert mu_out.mod["rna"].obsm["leiden_1"].shape[1] < mu_out.mod["rna"].obsp["distances"].shape[1], "Obsm should only contain a subset of the original obsp matrix"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
