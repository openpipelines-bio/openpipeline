import logging
import subprocess
from sys import stdout
import pytest
from pathlib import Path
from tempfile import NamedTemporaryFile
import anndata
import mudata
import numpy as np

## VIASH START
meta = {
    'executable': './target/docker/labels_transfer/knn/knn',
    'resources_dir': './resources_test/'
}
## VIASH END

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

@pytest.fixture
def reference_file():
    return f"{meta['resources_dir']}/annotation_test_data/TS_Blood_filtered.h5ad"

@pytest.fixture
def input_file():
    return f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"


def test_one_class(run_component, reference_file, input_file):
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
                "--query_obsm_key", "X_integrated_scanvi",
                "--reference", tempfile_reference_file.name,
                "--reference_obsm_key", "X_integrated_scanvi",
                "--output", "output.h5mu",
                "-k", "5",
                "--targets", "ann_level_1"])

    assert Path("output.h5mu").is_file()

    output_data = mudata.read_h5mu("output.h5mu")
    assert "ann_level_1_pred" in output_data.mod["rna"].obs
    assert "ann_level_1_uncertainty" in output_data.mod["rna"].obs
    assert "labels_transfer" in output_data.mod["rna"].uns
    assert "ann_level_1_pred" in output_data.mod["rna"].uns["labels_transfer"]
    assert output_data.mod["rna"].uns["labels_transfer"]["ann_level_1_pred"]["method"] == "KNN_pynndescent"
    assert output_data.mod["rna"].uns["labels_transfer"]["ann_level_1_pred"]["n_neighbors"] == 5
    assert "ann_level_2_pred" not in output_data.mod["rna"].obs

def test_multiple_classes(run_component, reference_file, input_file,):
    targets = ["ann_level_1", "ann_level_2", "ann_level_3", "ann_level_4", "ann_level_5", "ann_finest_level"]
    
    with NamedTemporaryFile("w", suffix=".h5ad") as tempfile_reference_file:
        reference_adata = anndata.read_h5ad(reference_file)
        reference_adata.obsm["X_integrated_scanvi"] = np.random.normal(size=(reference_adata.n_obs, 30))
        
        for i, target in enumerate(targets):
            class_names = [str(idx) for idx in range(i + 1)]  # e.g. ["0", "1", "2"], the higher the level, the more the classes
            reference_adata.obs[target] = np.random.choice(class_names, size=reference_adata.n_obs)

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
                "--query_obsm_key", "X_integrated_scanvi",
                "--reference", tempfile_reference_file.name,
                "--reference_obsm_key", "X_integrated_scanvi",
                "--output", "output.h5mu",
                "-k", "5",
                "--targets", ",".join(targets)])

    assert Path("output.h5mu").is_file()

    output_data = mudata.read_h5mu("output.h5mu")

    for target in targets:
        assert f"{target}_pred" in output_data.mod["rna"].obs
        assert f"{target}_uncertainty" in output_data.mod["rna"].obs
        assert "labels_transfer" in output_data.mod["rna"].uns
        assert f"{target}_pred" in output_data.mod["rna"].uns["labels_transfer"]
        assert output_data.mod["rna"].uns["labels_transfer"][f"{target}_pred"]["method"] == "KNN_pynndescent"
        assert output_data.mod["rna"].uns["labels_transfer"][f"{target}_pred"]["n_neighbors"] == 5
