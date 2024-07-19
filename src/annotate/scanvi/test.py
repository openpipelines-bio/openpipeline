import sys
import os
import pytest
import re
import mudata as mu
import anndata as ad
from openpipelinetestutils.asserters import assert_annotation_objects_equal
import scvi
import os
## VIASH START
meta = {
    "resources_dir": "resources_test"
}
## VIASH END

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"
reference_file = f"{meta['resources_dir']}/annotation_test_data/TS_Blood_filtered.h5ad"

@pytest.fixture
def create_scvi_model(random_path, tmp_path):
    def wrapper(input_file, reference_file):
        input_data = mu.read_h5mu(input_file)
        input_modality = input_data.mod["rna"]
        reference_data = ad.read_h5ad(reference_file)
        
        reference_data.var["gene_symbol"] = list(reference_data.var.index)
        reference_data.var.index = [re.sub("\\.[0-9]+$", "", s) for s in reference_data.var["ensemblid"]]

        common_ens_ids = list(set(reference_data.var.index).intersection(set(input_modality.var.index)))

        reference = reference_data[:, common_ens_ids].copy()
        query = input_modality[:, common_ens_ids].copy()

        scvi.model.SCVI.setup_anndata(reference,
                                    labels_key="cell_ontology_class"
                                    )

        scvi_model = scvi.model.SCVI(
            reference,
            use_layer_norm="both",
            use_batch_norm="none",
            encode_covariates=True,
            dropout_rate=0.2,
            n_layers=1,
            )
        scvi_model.train(max_epochs=10)
        
        input_data.mod["rna"] = query
        
        input_data_file = random_path(extension="h5mu")
        reference_file = random_path(extension="h5ad")
        scvi_model_file = tmp_path
        
        input_data.write_h5mu(input_data_file)
        reference.write_h5ad(reference_file)
        scvi_model.save(scvi_model_file, overwrite=True)
                
        return scvi_model_file, input_data_file, reference_file
    return wrapper

def test_simple_execution(run_component, random_h5mu_path, create_scvi_model):
    scvi_model_file, input_file_scvi, reference_file_scvi = create_scvi_model(input_file, reference_file)
    output_file = random_h5mu_path()

    run_component([
        "--input", input_file_scvi,
        "--reference", reference_file_scvi,
        "--scvi_reference_model", scvi_model_file,
        "--max_epochs", "10",
        "--output", output_file
    ])

    assert os.path.exists(output_file), "Output file does not exist"

    input_mudata = mu.read_h5mu(input_file_scvi)
    output_mudata = mu.read_h5mu(output_file)
    
    assert input_mudata.mod["rna"].n_obs == output_mudata.mod["rna"].n_obs, f"Number of observations changed"
    assert input_mudata.mod["rna"].n_vars == output_mudata.mod["rna"].n_vars, f"Number of variables changed"
    assert "X_scANVI" in output_mudata.mod["rna"].obsm, "Latent representation not added"
    assert "predictions_scanvi" in output_mudata.mod["rna"].obs, "Predictions not added"

    assert_annotation_objects_equal(input_mudata.mod["prot"],
                                    output_mudata.mod["prot"])
    
def test_multiple_arguments(run_component, random_h5mu_path, create_scvi_model, tmp_path):
    scvi_model_file, input_file_scvi, reference_file_scvi = create_scvi_model(input_file, reference_file)
    output_file = random_h5mu_path()

    run_component([
        "--input", input_file_scvi,
        "--reference", reference_file_scvi,
        "--scvi_reference_model", scvi_model_file,
        "--max_epochs", "10",
        "--output", output_file,
        "--input_obs_label", "scanvi_predictions",
        "--reduce_lr_on_plateau", "True",
        "--lr_patience", "5",
        "--lr_factor", "0.5",
        "--train_size", "0.8",
        "--early_stopping", "True",
        "--early_stopping_patience", "5",
        "--early_stopping_min_delta", "0.01",
        "--output_compression", "gzip",
        "--output_model", tmp_path
    ])

    assert os.path.exists(output_file), "Output file does not exist"
    assert os.path.exists(tmp_path / "model.pt"), "Model file does not exist" 

    input_mudata = mu.read_h5mu(input_file_scvi)
    output_mudata = mu.read_h5mu(output_file)
    
    assert input_mudata.mod["rna"].n_obs == output_mudata.mod["rna"].n_obs, f"Number of observations changed"
    assert input_mudata.mod["rna"].n_vars == output_mudata.mod["rna"].n_vars, f"Number of variables changed"
    assert "X_scANVI" in output_mudata.mod["rna"].obsm, "Latent representation not added"
    assert "scanvi_predictions" in output_mudata.mod["rna"].obs, "Predictions not added"

    assert_annotation_objects_equal(input_mudata.mod["prot"],
                                    output_mudata.mod["prot"])

if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))