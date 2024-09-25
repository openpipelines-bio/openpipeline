import sys
import os
import pytest
import subprocess
import re
import mudata as mu
from openpipelinetestutils.asserters import assert_annotation_objects_equal
import os
from sklearn.ensemble import RandomForestClassifier
import pickle

## VIASH START
meta = {
    "resources_dir": "resources_test"
}
## VIASH END

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"
reference_file = f"{meta['resources_dir']}/TS_Blood_filtered.h5mu"

@pytest.fixture
def subset_genes(random_h5mu_path):
    def wrapper(input_mudata_file, reference_mudata_file, modality):
        input_mudata = mu.read_h5mu(input_mudata_file)
        input_adata = input_mudata.mod[modality]
        reference_mudata = mu.read_h5mu(reference_mudata_file)
        reference_adata = reference_mudata.mod[modality]

        reference_mudata.var["gene_symbol"] = list(reference_mudata.var.index)
        reference_mudata.var.index = [re.sub("\\.[0-9]+$", "", s) for s in reference_mudata.var["ensemblid"]]
        reference_adata.var["gene_symbol"] = list(reference_adata.var.index)
        reference_adata.var.index = [re.sub("\\.[0-9]+$", "", s) for s in reference_adata.var["ensemblid"]]
        common_ens_ids = list(set(reference_adata.var.index).intersection(set(input_adata.var.index)))

        reference = reference_adata[:, common_ens_ids].copy()
        query = input_adata[:, common_ens_ids].copy()

        input_mudata.mod[modality] = query
        reference_mudata.mod[modality] = reference

        subset_input_mudata_file = random_h5mu_path()
        subset_reference_mudata_file = random_h5mu_path()

        input_mudata.write_h5mu(subset_input_mudata_file)
        reference_mudata.write_h5mu(subset_reference_mudata_file)
        return subset_input_mudata_file, subset_reference_mudata_file
    return wrapper

@pytest.fixture
def dummy_model(tmp_path, subset_genes):
    _, subset_reference_file = subset_genes(input_file, reference_file, "rna")
    reference_modality = mu.read_h5mu(subset_reference_file).mod["rna"]

    labels = reference_modality.obs["cell_ontology_class"].to_numpy()
    model = RandomForestClassifier()
    model.fit(reference_modality.X, labels)

    model_path = tmp_path / "model.pkl"
    with open(model_path, "wb") as f:
        pickle.dump(model, f)

    return model_path

def test_simple_execution(run_component, random_h5mu_path, subset_genes):
    subset_input_file, subset_reference_file = subset_genes(input_file, reference_file, "rna")
    output_file = random_h5mu_path()

    run_component([
        "--input", subset_input_file,
        "--reference", subset_reference_file,
        "--reference_obs_target", "cell_ontology_class",
        "--output", output_file
    ])

    assert os.path.exists(output_file), "Output file does not exist"

    input_mudata = mu.read_h5mu(input_file)
    output_mudata = mu.read_h5mu(output_file)

    assert_annotation_objects_equal(input_mudata.mod["prot"],
                                    output_mudata.mod["prot"])

    assert list(output_mudata.mod["rna"].obs.keys()) == ['random_forest_pred',
                                                         'random_forest_probability']

    obs_values = output_mudata.mod["rna"].obs["random_forest_probability"]
    assert all(0 <= value <= 1 for value in obs_values), "probabilities outside the range [0, 1]"

def test_custom_out_obs_model_params(run_component, random_h5mu_path, subset_genes):
    subset_input_file, subset_reference_file = subset_genes(input_file, reference_file, "rna")
    output_file = random_h5mu_path()

    run_component([
        "--input", subset_input_file,
        "--reference", subset_reference_file,
        "--reference_obs_target", "cell_ontology_class",
        "--output_obs_predictions", "dummy_pred",
        "--output_obs_probability", "dummy_probability",
        "--n_estimators", "10",
        "--criterion", "entropy",
        "--max_depth", "5",
        "--class_weight", "balanced",
        "--output", output_file
    ])

    assert os.path.exists(output_file), "Output file does not exist"

    input_mudata = mu.read_h5mu(input_file)
    output_mudata = mu.read_h5mu(output_file)

    assert_annotation_objects_equal(input_mudata.mod["prot"],
                                    output_mudata.mod["prot"])

    assert list(output_mudata.mod["rna"].obs.keys()) == ['dummy_pred',
                                                         'dummy_probability']

    obs_values = output_mudata.mod["rna"].obs["dummy_probability"]
    assert all(0 <= value <= 1 for value in obs_values), "probabilities outside the range [0, 1]"

def test_with_model(run_component, random_h5mu_path, dummy_model, subset_genes):
    subset_input_file, _ = subset_genes(input_file, reference_file, "rna")
    output_file = random_h5mu_path()

    run_component([
        "--input", subset_input_file,
        "--model", dummy_model,
        "--output", output_file,
        "--reference_obs_target", "cell_ontology_class"
    ])

    assert os.path.exists(output_file), "Output file does not exist"

    input_mudata = mu.read_h5mu(input_file)
    output_mudata = mu.read_h5mu(output_file)

    assert_annotation_objects_equal(input_mudata.mod["prot"],
                                    output_mudata.mod["prot"])

    assert list(output_mudata.mod["rna"].obs.keys()) == ['random_forest_pred',
                                                         'random_forest_probability']

    obs_values = output_mudata.mod["rna"].obs["random_forest_probability"]
    assert all(0 <= value <= 1 for value in obs_values), "probabilities outside the range [0, 1]"

def test_no_model_no_reference_error(run_component, random_h5mu_path):
    output_file = random_h5mu_path()

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component([
            "--input", input_file,
            "--output", output_file,
            "--reference_obs_target", "cell_ontology_class"
        ])
    assert re.search(r"ValueError: Make sure to provide either 'model' or 'reference', but not both.",
            err.value.stdout.decode('utf-8'))

def test_model_and_reference_error(run_component, random_h5mu_path, dummy_model, subset_genes):
    output_file = random_h5mu_path()
    subset_input_file, subset_reference_file = subset_genes(input_file, reference_file, "rna")
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component([
            "--input", subset_input_file,
            "--output", output_file,
            "--reference", subset_reference_file,
            "--reference_obs_target", "cell_ontology_class",
            "--model", dummy_model,
        ])
    assert re.search(r"ValueError: Make sure to provide either 'model' or 'reference', but not both.",
            err.value.stdout.decode('utf-8'))

def test_invalid_max_features(run_component, random_h5mu_path):
    output_file = random_h5mu_path()

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component([
            "--input", input_file,
            "--output", output_file,
            "--reference_obs_target", "cell_ontology_class",
            "--max_features", "invalid_value"
        ])
    assert re.search(r"Invaldid value invalid_value for --max_features: must either be an integer or one of 'sqrt', 'log2' or 'all'",
            err.value.stdout.decode('utf-8'))

if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))