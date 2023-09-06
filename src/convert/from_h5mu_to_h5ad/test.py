import sys
import pytest
import anndata as ad
import mudata as mu

## VIASH START
meta = {
    'executable': 'target/docker/convert/from_h5mu_to_h5ad/from_h5mu_to_h5ad',
    'resources_dir': 'resources_test'
}
## VIASH END

input = meta["resources_dir"] + "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"

def test_run(run_component, tmp_path):
    output = tmp_path / "output.h5ad"

    cmd_pars = [
        "--modality", "rna",
        "--input", input,
        "--output", str(output),
        "--output_compression", "gzip"
    ]
    run_component(cmd_pars)

    assert output.is_file(), "No output was created."

    adata = ad.read_h5ad(output)

    mdata = mu.read_h5mu(input)
    assert adata.n_obs == mdata.mod["rna"].n_obs
    assert adata.n_vars == mdata.mod["rna"].n_vars

if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))