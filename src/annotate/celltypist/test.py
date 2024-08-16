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

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"
reference_file = f"{meta['resources_dir']}/annotation_test_data/TS_Blood_filtered.h5mu"
model_file = f"{meta['resources_dir']}/annotation_test_data/celltypist_model_Immune_All_Low.pkl"
celltypist_input_file = f"{meta['resources_dir']}/annotation_test_data/demo_2000_cells.h5mu"

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
        "--reference_obs_targets", "cell_ontology_class",
        "--var_reference_gene_names", "ensemblid",
        "--output", output_file
    ])
    
    assert os.path.exists(output_file), "Output file does not exist"
    
    input_mudata = mu.read_h5mu(input_file_transformed)
    output_mudata = mu.read_h5mu(output_file)
    
    assert_annotation_objects_equal(input_mudata.mod["prot"],
                                    output_mudata.mod["prot"])
    
    assert {'celltypist_pred', 'celltypist_probability'}.issubset(output_mudata.mod["rna"].obs.keys()), "Required keys not found in .obs"
    
    obs_values = output_mudata.mod["rna"].obs["celltypist_probability"]
    assert all(0 <= value <= 1 for value in obs_values), ".obs at celltypist_probability has values outside the range [0, 1]"
    
def test_set_params(run_component, random_h5mu_path, normalize_log_transform):
    output_file = random_h5mu_path()
    input_file_transformed = normalize_log_transform(input_file, "rna")

    run_component([
        "--input", input_file_transformed,
        "--reference", reference_file,
        "--reference_obs_target", "cell_ontology_class",
        "--var_reference_gene_names", "ensemblid",
        "--feature_selection", "True",
        "--majority_voting", "True",
        "--C", "0.5",
        "--max_iter", "100",
        "--use_SGD",
        "--min_prop", "0.1",
        "--input_layer", "log_normalized",
        "--output", output_file,
        "--output_compression", "gzip",
    ])
    
    assert os.path.exists(output_file), "Output file does not exist"
    
    input_mudata = mu.read_h5mu(input_file_transformed)
    output_mudata = mu.read_h5mu(output_file)
    
    assert_annotation_objects_equal(input_mudata.mod["prot"],
                                    output_mudata.mod["prot"])
    
    assert {'celltypist_pred', 'celltypist_probability'}.issubset(output_mudata.mod["rna"].obs.keys()), "Required keys not found in .obs"
    
    obs_values = output_mudata.mod["rna"].obs["celltypist_probability"]
    assert all(0 <= value <= 1 for value in obs_values), ".obs at celltypist_probability has values outside the range [0, 1]"

def test_with_model(run_component, random_h5mu_path):
    output_file = random_h5mu_path()

    run_component([
        "--input", celltypist_input_file,
        "--model", model_file,
        "--reference_obs_targets", "cell_type",
        "--output", output_file
    ])
    
    assert os.path.exists(output_file), "Output file does not exist"
    
    output_mudata = mu.read_h5mu(output_file)
    
    assert {'celltypist_pred', 'celltypist_probability'}.issubset(output_mudata.mod["rna"].obs.keys()), "Required keys not found in .obs"
    
    obs_values = output_mudata.mod["rna"].obs["celltypist_probability"]
    assert all(0 <= value <= 1 for value in obs_values), ".obs at celltypist_probability has values outside the range [0, 1]"

def test_fail_check_reference_expression(run_component, random_h5mu_path):
    output_file = random_h5mu_path()

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component([
            "--input", input_file,
            "--reference", reference_file,
            "--var_reference_gene_names", "ensemblid",
            "--output", output_file,
            "--check_expression"
        ])
    assert re.search(r"Invalid expression matrix, expect log1p normalized expression to 10000 counts per cell",
            err.value.stdout.decode('utf-8'))
    
def test_fail_invalid_input_expression(run_component, random_h5mu_path):
    output_file = random_h5mu_path()

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component([
            "--input", input_file,
            "--reference", reference_file,
            "--var_reference_gene_names", "ensemblid",
            "--output", output_file
        ])
    assert re.search(r"Invalid expression matrix in `.X`, expect log1p normalized expression to 10000 counts per cell",
            err.value.stdout.decode('utf-8'))

if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))
        