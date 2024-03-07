import sys
import pytest
import subprocess
from mudata import read_h5mu
from openpipelinetestutils.asserters import assert_annotation_objects_equal

## VIASH START
meta = {
    'executable': './target/docker/dimred/tsne/tsne',
    'resources_dir': './resources_test/',
    'config': './src/dimred/tsne/config.vsh.yaml'
}
## VIASH END

input_path = meta["resources_dir"] + "pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"

def test_tsne(run_component, random_h5mu_path):
    output_path = random_h5mu_path()
    args = [
        "--input", input_path,
        "--output",  output_path,
        "--modality", "rna",
        "--use_rep", "X_pca",
        "--output_compression", "gzip"
    ]
    run_component(args)
    
    assert output_path.is_file(), "No output was created."
    output_mudata = read_h5mu(output_path)
    input_mudata = read_h5mu(input_path)
    
    # check whether tsne was found and remove for comparison
    assert "X_tsne" in output_mudata.mod["rna"].obsm, "Check whether output was found in .obsm"
    output_mudata.mod["rna"].obsm.pop("X_tsne")
    assert_annotation_objects_equal(output_mudata, input_mudata)

if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))