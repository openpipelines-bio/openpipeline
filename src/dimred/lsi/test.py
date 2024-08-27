import sys
import pytest
import subprocess
import mudata as mu
import numpy as np

## VIASH START
meta = {
    'resources_dir': 'resources_test',
    'executable': './target/docker/dimred/lsi/lsi',
    'config': './src/dimred/lsi/config.vsh.yaml'
}
## VIASH END

input_path = f"{meta['resources_dir']}/concat_test_data/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix_subset_unique_obs.h5mu"


'''
Tests: 
1. general test 
2. test HVF
3. test modality
4. test layer
5. test overwrite 
'''

@pytest.fixture
def atac_mudata(tmp_path):
    
    mdata = mu.read_h5mu(input_path)
    mdata.mod["atac"].layers["counts"] = mdata.mod["atac"].X
    mdata.mod["atac"].var["highly_variable"] = np.random.choice([True, False], size=mdata.mod["atac"].n_vars)
    print(mdata)

    mdata.write(tmp_path / "atac_mudata.h5mu")

    return tmp_path / "atac_mudata.h5mu"

# 1.general test
def test_lsi(tmp_path):
    output_path = tmp_path / "output_lsi.h5mu"
    # run component
    cmd_args = [
    meta["executable"],
     "--input", input_path,
     "--output", str(output_path),
     "--obsm_output", "X_test",
     "--num_components", "30"
    ]
    subprocess.run(cmd_args, check=True)
    

    assert output_path.is_file()
    data = mu.read_h5mu(output_path)
    assert "X_test" in data.mod['atac'].obsm
    assert data.mod["atac"].obsm["X_test"].shape == (data.mod["atac"].n_obs, 30)
    assert "lsi" in data.mod['atac'].uns
    assert "LSI" in data.mod['atac'].varm



# 2.test HVF 
def test_select_highly_variable_column(atac_mudata, tmp_path):
    output_path = tmp_path / "output_lsi.h5mu"

    # run component
    cmd_args = [
    meta["executable"],
     "--input", str(atac_mudata),
     "--output", str(output_path),
     "--var_input", "highly_variable"
    ]
    subprocess.run(cmd_args, check=True)
    
    assert output_path.is_file()
    data = mu.read_h5mu(output_path)
    assert "X_lsi" in data.mod['atac'].obsm
    assert data.mod["atac"].obsm["X_lsi"].shape == (data.mod["atac"].n_obs, 50)
    assert "highly_variable" in data.mod["atac"].var.columns
    assert "lsi" in data.mod['atac'].uns
    assert "LSI" in data.mod['atac'].varm
    assert data.mod["atac"].varm["LSI"].shape == (data.mod["atac"].n_vars, 50)


def test_highly_variable_column_does_not_exist_raises():
    with pytest.raises(subprocess.CalledProcessError) as err:
        cmd_args = [
        meta["executable"],
        "--input", input_path,
        "--output", "output_lsi.h5mu",
        "--var_input", "does_not_exist"
        ]
        subprocess.run(cmd_args, check=True)
        
        assert "ValueError: Requested to use .var column 'does_not_exist' as a selection of genes, but the column is not available." in \
            err.value.stdout.decode('utf-8')
        

# 3.test modality
def test_modality_does_not_exist_raises():
    with pytest.raises(subprocess.CalledProcessError) as err:
        cmd_args = [
        meta["executable"],
        "--input", input_path,
        "--output", "output_lsi.h5mu",
        "--modality", "does_not_exist"
        ]
        subprocess.run(cmd_args, check=True)
       
        assert "ValueError: Modality 'does_not_exist' was not found in mudata " + input_path +"." in \
            err.value.stdout.decode('utf-8')



# 4.test layer 
def test_selecting_input_layer(atac_mudata, tmp_path):
    output_path = tmp_path / "output_lsi.h5mu"

    # run component
    cmd_args = [
        meta["executable"],
        "--input", str(atac_mudata),
        "--output", str(output_path),
        "--num_components", "20",
        "--layer", "counts"
        ]
    subprocess.run(cmd_args, check=True)


    assert output_path.is_file()
    data = mu.read_h5mu(output_path)
    assert "counts" in data.mod["atac"].layers
    assert "X_lsi" in data.mod['atac'].obsm
    assert data.mod["atac"].obsm["X_lsi"].shape == (data.mod["atac"].n_obs, 20)
    assert "lsi" in data.mod['atac'].uns
    assert "LSI" in data.mod['atac'].varm



def test_raise_if_input_layer_is_missing():
    with pytest.raises(subprocess.CalledProcessError) as err:
        cmd_args = [
        meta["executable"],
        "--input", input_path,
        "--output", "output.h5mu",
        "--layer", "does_not_exist"
        ]
        subprocess.run(cmd_args, check=True)
        
        assert "ValueError: Layer 'does_not_exist' was not found in modality 'atac'." in \
            err.value.stdout.decode('utf-8')



# 5.test overwrite 

def test_output_field_already_present_raises(tmp_path):
    output_path = tmp_path / "output_lsi.h5mu"

    #create slots 
    input_data = mu.read_h5mu(input_path)
    input_data.mod["atac"].varm["LSI"] = np.zeros(shape=(input_data.mod["atac"].n_vars, 50))
    input_data.mod["atac"].obsm["X_lsi"] = np.zeros(shape=(input_data.mod["atac"].n_obs, 50))
    input_data.mod["atac"].uns['lsi'] = "test"
    tmp_file = tmp_path / "input_data_adjusted.h5mu"
    input_data.write_h5mu(tmp_file)

    with pytest.raises(subprocess.CalledProcessError) as err:
        cmd_args = [
        meta["executable"],
        "--input", str(tmp_file),
        "--output", str(output_path),
        "--output_compression", "gzip"
        ]
        subprocess.run(cmd_args, check=True)
       
        assert "ValueError: Requested to create field X_lsi in .obsm for " \
            "modality atac, but field already exists." in \
            err.value.stdout.decode('utf-8')

def test_output_field_already_present_overwrite(tmp_path):
    output_path = tmp_path / "output_lsi.h5mu"

    #create slots 
    input_data = mu.read_h5mu(input_path)
    input_data.mod["atac"].varm["LSI"] = np.zeros(shape=(input_data.mod["atac"].n_vars, 50))
    input_data.mod["atac"].obsm["X_lsi"] = np.zeros(shape=(input_data.mod["atac"].n_obs, 50))
    input_data.mod["atac"].uns['lsi'] = "test"
    tmp_file = tmp_path / "input_data_adjusted.h5mu"
    input_data.write_h5mu(tmp_file)

    cmd_args = [
        meta["executable"],
        "--input", str(tmp_file),
        "--output", str(output_path),
        "--output_compression", "gzip",
        "--overwrite",
         "--num_components", "30"
        ]
    subprocess.run(cmd_args, check=True)

    assert output_path.is_file()
    data = mu.read_h5mu(output_path)
    assert "X_lsi" in data.mod['atac'].obsm
    assert data.mod["atac"].obsm["X_lsi"].shape == (data.mod["atac"].n_obs, 30)
    assert "lsi" in data.mod['atac'].uns
    assert "LSI" in data.mod['atac'].varm

if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))