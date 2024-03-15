import sys
import pytest
import subprocess
from mudata import read_h5mu
from openpipelinetestutils.asserters import assert_annotation_objects_equal
import re

## VIASH START
meta = {
    'executable': './target/docker/dimred/densmap/densmap',
    'resources_dir': './resources_test/',
    'config': './src/dimred/densmap/config.vsh.yaml'
}
## VIASH END

input_path = meta["resources_dir"] + "pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"

# @pytest.fixture
# def mudata_no_obsm_pca(write_mudata_to_file):
#     input_mudata = read_h5mu(input_path)
#     input_mudata.mod["rna"].obsm.pop("X_pca")
#     return write_mudata_to_file(input_mudata)

def test_densmap(run_component, random_h5mu_path):
    output_path = random_h5mu_path()
    args = [
        "--input", input_path,
        "--output",  output_path,
        "--modality", "rna",
        "--output_compression", "gzip"
    ]
    run_component(args)
    
    assert output_path.is_file(), "No output was created."
    output_mudata = read_h5mu(output_path)
    input_mudata = read_h5mu(input_path)
    
    # check whether densmap was found and remove for comparison
    assert "X_densmap" in output_mudata.mod["rna"].obsm, "Check whether output was found in .obsm"
    assert "densmap" in output_mudata.mod["rna"].uns, "Check whether output was found in .uns"
    output_mudata.mod["rna"].obsm.pop("X_densmap")
    output_mudata.mod["rna"].uns.pop("densmap")
    assert_annotation_objects_equal(output_mudata, input_mudata)
    
def test_densmap_custom_obsm_output(run_component, random_h5mu_path):
    output_path = random_h5mu_path()
    args = [
        "--input", input_path,
        "--output",  output_path,
        "--modality", "rna",
        "--output_compression", "gzip",
        "--obsm_output", "X_custom_densmap"
    ]
    run_component(args)
    
    assert output_path.is_file(), "No output was created."
    output_mudata = read_h5mu(output_path)
    input_mudata = read_h5mu(input_path)
    
    # check whether tsne was found and remove for comparison
    assert "X_custom_densmap" in output_mudata.mod["rna"].obsm, "Check whether output was found in .obsm"
    assert "densmap" in output_mudata.mod["rna"].uns, "Check whether output was found in .uns"
    output_mudata.mod["rna"].obsm.pop("X_custom_densmap")
    output_mudata.mod["rna"].uns.pop("densmap")
    assert_annotation_objects_equal(output_mudata, input_mudata)
    
def test_densmap_no_neighbors_raise(run_component, random_h5mu_path):
    output_path = random_h5mu_path()
    args = [
        "--input", input_path,
        "--output",  output_path,
        "--modality", "prot",
        "--output_compression", "gzip"
    ]
    
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(args)
    assert re.search(r"ValueError: 'neighbors' was not found in .mod\['prot'\].uns.",
                     err.value.stdout.decode('utf-8'))
        
    
if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))