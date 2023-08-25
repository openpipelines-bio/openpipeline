import pytest
from pathlib import Path
from tempfile import NamedTemporaryFile
import anndata
import mudata
import numpy as np

## VIASH START
meta = {
    'executable': './target/docker/labels_transfer/xgboost/xgboost',
    'resources_dir': './resources_test/'
}
## VIASH END

@pytest.fixture
def reference_file():
    return f"{meta['resources_dir']}/annotation_test_data/TS_Blood_filtered.h5ad"

@pytest.fixture
def input_file():
    return f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"

def test_one_class(run_component, input_file, reference_file):
    with NamedTemporaryFile("w", suffix=".h5ad") as tempfile_reference_file:
        reference_adata = anndata.read_h5ad(reference_file)
        reference_adata.obsm["X_integrated_scanvi"] = np.random.normal(size=(reference_adata.n_obs, 30))
        reference_adata.obs["ann_level_1"] = np.random.choice(["CD4 T cells", "CD8 T cells"], size=reference_adata.n_obs)
        reference_adata.write(tempfile_reference_file.name)
        
        with NamedTemporaryFile("w", suffix=".h5mu") as tempfile_input_file:
            input_data = mudata.read_h5mu(input_file)
            adata = input_data.mod["rna"]

            # Simulate a latent embedding
            # Later, we can use a real one by calling `integrate` component in advamnce
            adata.obsm["X_integrated_scanvi"] = np.random.normal(size=(adata.n_obs, 30))

            input_data.write(tempfile_input_file.name)

            run_component([
                "--input", tempfile_input_file.name,
                "--modality", "rna",
                "--input_obsm_features", "X_integrated_scanvi",
                "--reference", tempfile_reference_file.name,
                "--reference_obsm_features", "X_integrated_scanvi",
                "--reference_obs_targets", "ann_level_1",
                "--output", "output.h5mu",
                "--model_output", "model_one_class",
                "--use_gpu", "false",
                "--max_depth", "6"
            ])

    assert Path("output.h5mu").is_file()

    output_data = mudata.read_h5mu("output.h5mu")
    assert "ann_level_1_pred" in output_data.mod["rna"].obs, f"Predictions are missing from output\noutput: {output_data.mod['rna'].obs}"
    assert "ann_level_1_uncertainty" in output_data.mod["rna"].obs, f"Uncertainties are missing from output\noutput: {output_data.mod['rna'].obs}"
    assert "labels_transfer" in output_data.mod["rna"].uns, f"Parameters are missing from output\noutput: {output_data.mod['rna'].uns}"
    assert "ann_level_1" in output_data.mod["rna"].uns["labels_transfer"], f"Parameters are missing from output\noutput: {output_data.mod['rna'].uns}"
    assert output_data.mod["rna"].uns["labels_transfer"]["ann_level_1"].get("method") == "XGBClassifier", f"Wrong method in parameters\noutput: {output_data.mod['rna'].uns}"
    assert output_data.mod["rna"].uns["labels_transfer"]["ann_level_1"].get("max_depth") == 6, f"Wrong max depth in parameters\noutput: {output_data.mod['rna'].uns}"
    assert "ann_level_2" not in output_data.mod["rna"].obs, f"Level 2 should not be in output\noutput: {output_data.mod['rna'].obs}"

    # Remove output file to prevent errors in the next tests
    Path("output.h5mu").unlink()

def test_multiple_classes(run_component, input_file, reference_file):
    targets = ["ann_level_1", "ann_level_2", "ann_level_3", "ann_level_4", "ann_level_5", "ann_finest_level"]

    with NamedTemporaryFile("w", suffix=".h5ad") as tempfile_reference_file:
        reference_adata = anndata.read_h5ad(reference_file)
        reference_adata.obsm["X_integrated_scvi"] = np.random.normal(size=(reference_adata.n_obs, 30))
        for i, target in enumerate(targets):
            class_names = [str(idx) for idx in range(i + 2)]  # e.g. ["0", "1", "2"], the higher the level, the more the classes
            reference_adata.obs[target] = np.random.choice(class_names, size=reference_adata.n_obs)

        reference_adata.write(tempfile_reference_file.name)

        with NamedTemporaryFile("w", suffix=".h5mu") as tempfile_input_file:
            input_data = mudata.read_h5mu(input_file)
            adata = input_data.mod["rna"]

            # Simulate a latent embedding
            # Later, we can use a real one by calling `integrate` component in advamnce
            adata.obsm["X_integrated_scvi"] = np.random.normal(size=(adata.n_obs, 30))

            input_data.write(tempfile_input_file.name)

            run_component([
                "--input", tempfile_input_file.name,
                "--modality", "rna",
                "--input_obsm_features", "X_integrated_scvi",
                "--reference", tempfile_reference_file.name,
                "--reference_obsm_features", "X_integrated_scvi",
                "--reference_obs_targets", ",".join(targets),
                "--output", "output.h5mu",
                "--model_output", "model_multiple_classes",
                "--max_depth", "6",
                "--use_gpu", "false"
            ])

    assert Path("output.h5mu").is_file()

    output_data = mudata.read_h5mu("output.h5mu")

    for target in targets:
        assert f"{target}_pred" in output_data.mod["rna"].obs, f"Predictions are missing from output\noutput: {output_data.mod['rna'].obs}"
        assert f"{target}_uncertainty" in output_data.mod["rna"].obs, f"Uncertainties are missing from output\noutput: {output_data.mod['rna'].obs}"
        assert "labels_transfer" in output_data.mod["rna"].uns, f"Parameters are missing from output\noutput: {output_data.mod['rna'].uns}"
        assert target in output_data.mod["rna"].uns["labels_transfer"], f"Parameters are missing from output\noutput: {output_data.mod['rna'].uns}"
        assert output_data.mod["rna"].uns["labels_transfer"][target].get("method") == "XGBClassifier", f"Wrong method in parameters\noutput: {output_data.mod['rna'].uns}"
        assert output_data.mod["rna"].uns["labels_transfer"][target].get("max_depth") == 6, f"Wrong max depth in parameters\noutput: {output_data.mod['rna'].uns}"

    # Remove output file to prevent errors in the next tests
    Path("output.h5mu").unlink()

def test_retraining(run_component, input_file, reference_file):
    with NamedTemporaryFile("w", suffix=".h5ad") as tempfile_reference_file:
        reference_adata = anndata.read_h5ad(reference_file)
        reference_adata.obsm["X_integrated_scanvi"] = np.random.normal(size=(reference_adata.n_obs, 30))
        reference_adata.obs["ann_level_1"] = np.random.choice(["CD4 T cells", "CD8 T cells"], size=reference_adata.n_obs)
        reference_adata.obs["ann_level_2"] = np.random.choice(["CD4 naive T cells", "CD4 memory T cells", "CD8 naive T cells", "CD8 memory T cells"], size=reference_adata.n_obs)
        reference_adata.obs["ann_level_3"] = np.random.choice(["CD4 naive T cells", "CD4 memory T cells", "CD8 naive T cells", "CD8 memory T cells", "T Rex cells"], size=reference_adata.n_obs)
        reference_adata.write(tempfile_reference_file.name)

        targets = ["ann_level_1", "ann_level_2"]
    
        with NamedTemporaryFile("w", suffix=".h5mu") as tempfile_input_file:
            input_data = mudata.read_h5mu(input_file)
            adata = input_data.mod["rna"]

            # Simulate a latent embedding
            # Later, we can use a real one by calling `integrate` component in advance
            adata.obsm["X_integrated_scanvi"] = np.random.normal(size=(adata.n_obs, 30))

            input_data.write(tempfile_input_file.name)

            # Train first 2 targets
            run_component([
                "--input", tempfile_input_file.name,
                "--modality", "rna",
                "--input_obsm_features", "X_integrated_scanvi",
                "--reference", tempfile_reference_file.name,
                "--reference_obsm_features", "X_integrated_scanvi",
                "--output", "output.h5mu",
                "--model_output", "model_retraining",
                "--max_depth", "6",
                "--reference_obs_targets", ",".join(targets)])

            assert Path("output.h5mu").is_file()

            output_data = mudata.read_h5mu("output.h5mu")

            for target in targets:
                assert f"{target}_pred" in output_data.mod["rna"].obs, f"Predictions are missing from output\noutput: {output_data.mod['rna'].obs}"
                assert f"{target}_uncertainty" in output_data.mod["rna"].obs, f"Uncertainties are missing from output\noutput: {output_data.mod['rna'].obs}"
                assert "labels_transfer" in output_data.mod["rna"].uns, f"Parameters are missing from output\noutput: {output_data.mod['rna'].uns}"
                assert target in output_data.mod["rna"].uns["labels_transfer"], f"Parameters are missing from output\noutput: {output_data.mod['rna'].uns}"
                assert output_data.mod["rna"].uns["labels_transfer"][target].get("method") == "XGBClassifier", f"Wrong method in parameters\noutput: {output_data.mod['rna'].uns}"
                assert output_data.mod["rna"].uns["labels_transfer"][target].get("max_depth") == 6, f"Wrong max depth in parameters\noutput: {output_data.mod['rna'].uns}"

            # Add one more target
            # Now the code should use 2 previously trained models, and train only the third
            targets.append("ann_level_3")
            run_component([
                "--input", "output.h5mu",
                "--modality", "rna",
                "--input_obsm_features", "X_integrated_scanvi",
                "--reference", tempfile_reference_file.name,
                "--reference_obsm_features", "X_integrated_scanvi",
                "--model_output", "model_retraining",
                "--output", "output.h5mu",
                "--max_depth", "4",  # Change parameter so that we could make sure if the model is old or new
                "--reference_obs_targets", ",".join(targets)])

    assert Path("output.h5mu").is_file()

    output_data = mudata.read_h5mu("output.h5mu")

    for target in targets:
        assert f"{target}_pred" in output_data.mod["rna"].obs, f"Predictions are missing from output\noutput: {output_data.mod['rna'].obs}"
        assert f"{target}_uncertainty" in output_data.mod["rna"].obs, f"Uncertainties are missing from output\noutput: {output_data.mod['rna'].obs}"
        assert "labels_transfer" in output_data.mod["rna"].uns, f"Parameters are missing from output\noutput: {output_data.mod['rna'].uns}"
        assert target in output_data.mod["rna"].uns["labels_transfer"]
        assert output_data.mod["rna"].uns["labels_transfer"][target].get("method") == "XGBClassifier"
    
    assert output_data.mod["rna"].uns["labels_transfer"][targets[0]].get("max_depth") == 6
    assert output_data.mod["rna"].uns["labels_transfer"][targets[1]].get("max_depth") == 6
    assert output_data.mod["rna"].uns["labels_transfer"][targets[2]].get("max_depth") == 4

    # Remove output file to prevent errors in the next tests
    Path("output.h5mu").unlink()


if __name__ == '__main__':
    exit(pytest.main([__file__]))