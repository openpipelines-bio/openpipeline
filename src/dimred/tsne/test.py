import sys
import pytest
import subprocess
from mudata import read_h5mu
from openpipelinetestutils.asserters import assert_annotation_objects_equal
import re

## VIASH START
meta = {
    'executable': './target/docker/dimred/tsne/tsne',
    'resources_dir': './resources_test/',
    'config': './src/dimred/tsne/config.vsh.yaml'
}
## VIASH END

input_path = meta["resources_dir"] + "pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"

@pytest.fixture
def mudata_no_obsm_pca(write_mudata_to_file):
    input_mudata = read_h5mu(input_path)
    input_mudata.mod["rna"].obsm.pop("X_pca")
    return write_mudata_to_file(input_mudata)

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
    

def test_tsne_no_pca_in_input(run_component, random_h5mu_path, mudata_no_obsm_pca):
    output_path = random_h5mu_path()
    run_component([
        "--input", mudata_no_obsm_pca,
        "--output",  output_path,
        "--modality", "rna",
        "--n_pcs", "100",
        "--output_compression", "gzip"
    ])
    
    assert output_path.is_file(), "No output was created."
    output_mudata = read_h5mu(output_path)
    input_mudata = read_h5mu(mudata_no_obsm_pca)
    
    # check whether tsne was found and remove for comparison
    assert "X_tsne" in output_mudata.mod["rna"].obsm, "Check whether output was found in .obsm"
    output_mudata.mod["rna"].obsm.pop("X_tsne")
    assert_annotation_objects_equal(output_mudata, input_mudata)
    
    
def test_tsne_too_many_pcs_raise(run_component, random_h5mu_path):
    output_path = random_h5mu_path()
    args = [
        "--input", input_path,
        "--output",  output_path,
        "--modality", "rna",
        "--use_rep", "X_pca",
        "--output_compression", "gzip",
        "--n_pcs", "100"
    ]
    
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(args)
    assert re.search(r"ValueError: X_pca does not have enough Dimensions\. Provide a Representation with equal or more dimensions than`n_pcs` or lower `n_pcs`",
                     err.value.stdout.decode('utf-8'))
    
    
def test_tsne_raw_sparse_X_raise(run_component, random_h5mu_path):
    output_path = random_h5mu_path()
    args = [
        "--input", input_path,
        "--output",  output_path,
        "--modality", "rna",
        "--use_rep", "X",
        "--output_compression", "gzip"
    ]
    
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(args)
    assert re.search(r"TypeError: PCA initialization is currently not supported with the sparse input matrix\.",
                     err.value.stdout.decode('utf-8'))
    
    
def test_tsne_no_pca_in_input_raise(run_component, random_h5mu_path, mudata_no_obsm_pca):
    output_path = random_h5mu_path()
    args = [
            "--input", mudata_no_obsm_pca,
            "--output",  output_path,
            "--modality", "rna",
            "--use_rep", "X_pca",
            "--output_compression", "gzip"
    ]
    
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(args)
    assert re.search(r"ValueError: Did not find X_pca in `\.obsm.keys\(\)`\. You need to compute it first\.",
                     err.value.stdout.decode('utf-8'))
    

if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))