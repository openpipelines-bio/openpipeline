import pytest
from pathlib import Path
from tempfile import NamedTemporaryFile

import mudata
from anndata.tests.helpers import assert_equal

## VIASH START
meta = {
    'executable': './target/docker/integrate/scvi/scvi',
    'resources_dir': './resources_test/pbmc_1k_protein_v3/'
}
## VIASH END

import sys
sys.path.append(meta['resources_dir'])
from subset_vars import subset_vars

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3_mms.h5mu"

def test_scvi(run_component):
    with NamedTemporaryFile('w', suffix=".h5mu") as tempfile_input_file:
        input_data = mudata.read_h5mu(input_file)
        input_rna = input_data.mod['rna']
        input_rna.obs['batch'] = 'A'
        column_index = input_rna.obs.columns.get_indexer(['batch'])
        input_rna.obs.iloc[slice(input_rna.n_obs//2, None), column_index] = 'B'
        input_data.write(tempfile_input_file.name)

        run_component([
            "--input", tempfile_input_file.name,
            "--modality", "rna",
            "--obs_batch", "batch",
            "--output", "output.h5mu",
            "--model_output", "test/",
            "--max_epochs", "1",
            "--n_obs_min_count", "10",
            "--n_var_min_count", "10",
            "--output_compression", "gzip"])
        
        # check files
        assert Path("output.h5mu").is_file(), "Output file does not exist"
        assert Path("test").is_dir()
        assert Path("test/model.pt").is_file()

        # check output h5mu
        output_data = mudata.read_h5mu("output.h5mu")
        output_rna = output_data.mod['rna']
        assert output_rna.n_obs == input_rna.n_obs, f"Number of observations changed\noutput_data: {output_data}"
        assert output_rna.n_vars == input_rna.n_vars, f"Number of variables changed\noutput_data: {output_data}"
        assert 'X_scvi_integrated' in output_rna.obsm, f"X_scvi_integrated not added\noutput_data: {output_data}"

        # assert that nothing else has changed
        del output_rna.obsm["X_scvi_integrated"]
        assert_equal(input_rna, output_rna)

def test_scvi_after_hvg(run_component):
    with NamedTemporaryFile('w', suffix=".h5mu") as tempfile_input_file:
        input_data = mudata.read_h5mu(input_file)
        input_rna = input_data.mod['rna']
        input_rna.obs['batch2'] = 'A'
        column_index = input_rna.obs.columns.get_indexer(['batch2'])
        input_rna.obs.iloc[slice(input_rna.n_obs//2, None), column_index] = 'B'
        input_data.write(tempfile_input_file.name)

        run_component([
            "--input", tempfile_input_file.name,
            "--modality", "rna",
            "--obs_batch", "batch2",
            "--output", "output.h5mu",
            "--model_output", "test/",
            "--var_input", "filter_with_hvg",
            "--max_epochs", "1",
            "--n_obs_min_count", "10",
            "--n_var_min_count", "10",
            "--output_compression", "gzip",
            "--obsm_output", "X_scvi_integrated2"])
        
        # check files
        assert Path("output.h5mu").is_file(), "Output file does not exist"
        assert Path("test").is_dir()
        assert Path("test/model.pt").is_file()

        # check output h5mu
        output_data = mudata.read_h5mu("output.h5mu")
        output_rna = output_data.mod['rna']
        assert output_rna.n_obs == input_rna.n_obs, f"Number of observations changed\noutput_data: {output_data}"
        assert output_rna.n_vars == input_rna.n_vars, f"Number of variables changed\noutput_data: {output_data}"
        assert 'X_scvi_integrated2' in output_rna.obsm, f"X_scvi_integrated not added\noutput_data: {output_data}"

        # assert that nothing else has changed
        del output_rna.obsm["X_scvi_integrated2"]
        assert_equal(input_rna, output_rna)

def test_hvg_subsetting_helper():
    input_data = mudata.read_h5mu(input_file)
    adata = input_data.mod["rna"]

    old_n_genes = adata.n_vars

    adata.var["highly_variable_features"] = False
    adata.var.iloc[:old_n_genes // 2, adata.var.columns.get_indexer(["highly_variable_features"])] = True

    adata = subset_vars(adata, subset_col="highly_variable_features")

    # Correct number of genes is subsetted
    assert adata.n_vars == old_n_genes // 2
    # Only HVG are subsetted
    assert adata.var["highly_variable_features"].all()
       
if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))