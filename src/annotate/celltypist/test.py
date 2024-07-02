import sys
import os
import pytest
import subprocess
import re
import mudata as mu
import scanpy as sc
import anndata as ad
from openpipelinetestutils.asserters import assert_annotation_objects_equal
## VIASH START
meta = {
    "resources_dir": "resources_test"
}
## VIASH END

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"
reference_file = f"{meta['resources_dir']}/annotation_test_data/TS_Blood_filtered.h5ad"
model_file = f"{meta['resources_dir']}/annotation_test_data/celltypist_model_blood_filtered.pkl"

@pytest.fixture
def normalize_log_transform(random_h5mu_path):
    def wrapper(input_mudata_file, modality, target_sum=1e4):
        input_mudata = mu.read_h5mu(input_mudata_file)
        input_adata = input_mudata.mod[modality]
        adata = input_adata.copy()
        input_layer = adata.X
        data_for_scanpy = ad.AnnData(X=input_layer.copy())
        sc.pp.normalize_total(data_for_scanpy, target_sum=target_sum)
        sc.pp.log1p(data_for_scanpy,
                    base=None,
                    layer=None, # use X
                    copy=False) # allow overwrites in the copy that was made
        adata.X = data_for_scanpy.X
        adata.uns['log1p'] = data_for_scanpy.uns['log1p'].copy()
        input_mudata.mod[modality] = adata
        transformed_input_mudata_file = random_h5mu_path()
        input_mudata.write_h5mu(transformed_input_mudata_file)
        return transformed_input_mudata_file
    return wrapper

def test_simple_execution(run_component, random_h5mu_path, normalize_log_transform):
    output_file = random_h5mu_path()
    input_file_transformed = normalize_log_transform(input_file, "rna")

    run_component([
        "--input", input_file_transformed,
        "--reference", reference_file,
        "--output", output_file
    ])
    
    assert os.path.exists(output_file), "Output file does not exist"
    
    input_mudata = mu.read_h5mu(input_file_transformed)
    output_mudata = mu.read_h5mu(output_file)
    
    assert_annotation_objects_equal(input_mudata.mod["prot"],
                                    output_mudata.mod["prot"])
    
    assert list(output_mudata.mod["rna"].obs.keys()) == ['predicted_labels',
                                                         'over_clustering',
                                                         'majority_voting',
                                                         'conf_score']

def test_with_model(run_component, random_h5mu_path, normalize_log_transform):
    output_file = random_h5mu_path()
    input_file_transformed = normalize_log_transform(input_file, "rna")

    run_component([
        "--input", input_file_transformed,
        "--model", model_file,
        "--output", output_file
    ])
    
    assert os.path.exists(output_file), "Output file does not exist"
    
    input_mudata = mu.read_h5mu(input_file_transformed)
    output_mudata = mu.read_h5mu(output_file)
    
    assert_annotation_objects_equal(input_mudata.mod["prot"],
                                    output_mudata.mod["prot"])
    
    assert list(output_mudata.mod["rna"].obs.keys()) == ['predicted_labels',
                                                         'over_clustering',
                                                         'majority_voting',
                                                         'conf_score']

def test_more_params_save_model(run_component, random_h5mu_path, tmp_path, normalize_log_transform):
    output_file = random_h5mu_path()
    input_file_transformed = normalize_log_transform(input_file, "rna")
    model_save_path = tmp_path / "new_model.pkl"

    run_component([
        "--input", input_file_transformed,
        "--reference", reference_file,
        "--modality", "rna",
        "--reference_obs_label", "cell_ontology_class",
        "--output", output_file,
        "--majority_voting", "False",
        "--model_save_path", model_save_path
    ])
    
    assert os.path.exists(output_file), "Output file does not exist"
    assert os.path.exists(model_save_path), "Model file does not exist"
    
    input_mudata = mu.read_h5mu(input_file_transformed)
    output_mudata = mu.read_h5mu(output_file)
    
    assert_annotation_objects_equal(input_mudata.mod["prot"],
                                    output_mudata.mod["prot"])
    
    assert list(output_mudata.mod["rna"].obs.keys()) == ['predicted_labels',
                                                         'conf_score']

def test_fail_check_reference_expression(run_component, random_h5mu_path, normalize_log_transform):
    output_file = random_h5mu_path()

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component([
            "--input", input_file,
            "--reference", reference_file,
            "--output", output_file,
            "--check_expression", "True"
        ])
    assert re.search(r"ValueError: 🛑 Invalid expression matrix, expect log1p normalized expression to 10000 counts per cell",
            err.value.stdout.decode('utf-8'))
    
def test_fail_invalid_input_expression(run_component, random_h5mu_path):
    output_file = random_h5mu_path()

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component([
            "--input", input_file,
            "--reference", reference_file,
            "--output", output_file
        ])
    assert re.search(r"ValueError: 🛑 Invalid expression matrix in `.X`, expect log1p normalized expression to 10000 counts per cell",
            err.value.stdout.decode('utf-8'))

if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))
        