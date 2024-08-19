import sys
import pytest
import mudata
import numpy as np

## VIASH START
meta = {
    'executable': './target/executable/interpret/lianapy/',
    'resources_dir': './resources_test/pbmc_1k_protein_v3/'
}
## VIASH END

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3_mms.h5mu"


def test_lianapy(run_component, tmp_path):
    output_path = tmp_path / "output.h5mu"

    # run component
    run_component([
        "--input", input_file,
        "--output_compression", "gzip",
        "--modality", "rna",
        "--layer", "log_normalized",
        "--groupby", "harmony_integration_leiden_1.0",
        "--resource_name", "consensus",
        "--gene_symbol", "gene_symbol",
        "--expr_prop", "0.1",
        "--min_cells", "5",
        "--aggregate_method", "rra",
        "--return_all_lrs", "False",
        "--n_perms", "11",
        "--output", str(output_path)])
    assert output_path.is_file()

    # check output
    input_data = mudata.read_h5mu(input_file)
    output_data = mudata.read_h5mu(output_path)
    np.testing.assert_array_equal(output_data.mod['rna'].X.data, input_data.mod['rna'].X.data)
    np.testing.assert_array_equal(input_data.mod['rna'].var.index, output_data.mod['rna'].var.index)
    assert "liana_res" in output_data.mod["rna"].uns
    assert all(elem in output_data.mod['rna'].obs['harmony_integration_leiden_1.0'].values for elem in output_data.mod['rna'].uns['liana_res']['source'].unique())
    assert all(elem in output_data.mod['rna'].obs['harmony_integration_leiden_1.0'].values for elem in output_data.mod['rna'].uns['liana_res']['target'].unique())

if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))