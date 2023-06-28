import unittest
import subprocess
from pathlib import Path
from tempfile import NamedTemporaryFile
import os

import anndata
import mudata
import numpy as np


## VIASH START
meta = {
    'executable': './target/docker/labels_transfer/xgboost/xgboost',
    'resources_dir': './resources_test/'
}
## VIASH END

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"
reference_file = f"{meta['resources_dir']}/annotation_test_data/TS_Blood_filtered.h5ad"

class TestXGBoost(unittest.TestCase):
    def _run_and_check_output(self, args_as_list):
        try:
            subprocess.check_output([meta['executable']] + args_as_list)
        except subprocess.CalledProcessError as e:
            print(e.stdout.decode("utf-8"))
            raise e

    def test_one_class(self):

        with NamedTemporaryFile("w", suffix=".h5ad") as tempfile_reference_file:
            reference_adata = anndata.read_h5ad(reference_file)
            reference_adata.obsm["X_integrated_scvi"] = np.random.normal(size=(reference_adata.n_obs, 30))
            reference_adata.obs["ann_level_1"] = np.random.choice(["CD4 T cells", "CD8 T cells"], size=reference_adata.n_obs)
            reference_adata.write(tempfile_reference_file.name)
        
            with NamedTemporaryFile("w", suffix=".h5mu") as tempfile_input_file:
                input_data = mudata.read_h5mu(input_file)
                adata = input_data.mod["rna"]

                # Simulate a latent embedding
                # Later, we can use a real one by calling `integrate` component in advamnce
                adata.obsm["X_integrated_scvi"] = np.random.normal(size=(adata.n_obs, 30))

                input_data.write(tempfile_input_file.name)

                self._run_and_check_output([
                    "--input", tempfile_input_file.name,
                    "--modality", "rna",
                    "--query_obsm_key", "X_integrated_scvi",
                    "--reference", tempfile_reference_file.name,
                    "--reference_obsm_key", "X_integrated_scvi",
                    "--model_output", "model_one_class",
                    "--output", "output.h5mu",
                    "--targets", "ann_level_1"])

        self.assertTrue(Path("output.h5mu").is_file())

        output_data = mudata.read_h5mu("output.h5mu")
        self.assertIn("ann_level_1_pred", output_data.mod["rna"].obs)
        self.assertIn("ann_level_1_pred_uncertainty", output_data.mod["rna"].obs)
        self.assertIn("labels_transfer", output_data.mod["rna"].uns)
        self.assertIn("ann_level_1_pred", output_data.mod["rna"].uns["labels_transfer"])
        self.assertEqual("XGBClassifier", output_data.mod["rna"].uns["labels_transfer"]["ann_level_1_pred"]["method"])
        self.assertEqual(6, output_data.mod["rna"].uns["labels_transfer"]["ann_level_1_pred"]["max_depth"])
        self.assertNotIn("ann_level_2_pred", output_data.mod["rna"].obs)

        # Remove output file to prevent errors in the next tests
        Path("output.h5mu").unlink()

    def test_multiple_classes(self):

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

                self._run_and_check_output([
                    "--input", tempfile_input_file.name,
                    "--modality", "rna",
                    "--query_obsm_key", "X_integrated_scanvi",
                    "--reference", tempfile_reference_file.name,
                    "--reference_obsm_key", "X_integrated_scanvi",
                    "--output", "output.h5mu",
                    "--model_output", "model_multiple_classes",
                    "--max_depth", "6",
                    "--targets", ",".join(targets)])

        self.assertTrue(Path("output.h5mu").is_file())

        output_data = mudata.read_h5mu("output.h5mu")

        for target in targets:
            self.assertIn(f"{target}_pred", output_data.mod["rna"].obs)
            self.assertIn(f"{target}_pred_uncertainty", output_data.mod["rna"].obs)
            self.assertIn("labels_transfer", output_data.mod["rna"].uns)
            self.assertIn(f"{target}_pred", output_data.mod["rna"].uns["labels_transfer"])
            self.assertEqual("XGBClassifier", output_data.mod["rna"].uns["labels_transfer"][f"{target}_pred"]["method"])
            self.assertEqual(6, output_data.mod["rna"].uns["labels_transfer"][f"{target}_pred"]["max_depth"])

        # Remove output file to prevent errors in the next tests
        Path("output.h5mu").unlink()

    def test_retraining(self):

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
                self._run_and_check_output([
                    "--input", tempfile_input_file.name,
                    "--modality", "rna",
                    "--query_obsm_key", "X_integrated_scanvi",
                    "--reference", tempfile_reference_file.name,
                    "--reference_obsm_key", "X_integrated_scanvi",
                    "--output", "output.h5mu",
                    "--model_output", "model_retraining",
                    "--max_depth", "6",
                    "--targets", ",".join(targets)])

                self.assertTrue(Path("output.h5mu").is_file())

                output_data = mudata.read_h5mu("output.h5mu")

                for target in targets:
                    self.assertIn(f"{target}_pred", output_data.mod["rna"].obs)
                    self.assertIn(f"{target}_pred_uncertainty", output_data.mod["rna"].obs)
                    self.assertIn("labels_transfer", output_data.mod["rna"].uns)
                    self.assertIn(f"{target}_pred", output_data.mod["rna"].uns["labels_transfer"])
                    self.assertEqual("XGBClassifier", output_data.mod["rna"].uns["labels_transfer"][f"{target}_pred"]["method"])
                    self.assertEqual(6, output_data.mod["rna"].uns["labels_transfer"][f"{target}_pred"]["max_depth"])

                # Add one more target
                # Now the code should use 2 previously trained models, and train only the third
                targets.append("ann_level_3")
                self._run_and_check_output([
                    "--input", tempfile_input_file.name,
                    "--modality", "rna",
                    "--query_obsm_key", "X_integrated_scanvi",
                    "--reference", tempfile_reference_file.name,
                    "--reference_obsm_key", "X_integrated_scanvi",
                    "--model_output", "model_retraining",
                    "--output", "output.h5mu",
                    "--max_depth", "4",  # Change parameter so that we could make sure if the model is old or new
                    "--targets", ",".join(targets)])

        self.assertTrue(Path("output.h5mu").is_file())

        output_data = mudata.read_h5mu("output.h5mu")

        for target in targets:
            self.assertIn(f"{target}_pred", output_data.mod["rna"].obs)
            self.assertIn(f"{target}_pred_uncertainty", output_data.mod["rna"].obs)
            self.assertIn("labels_transfer", output_data.mod["rna"].uns)
            self.assertIn(f"{target}_pred", output_data.mod["rna"].uns["labels_transfer"])
            self.assertEqual("XGBClassifier", output_data.mod["rna"].uns["labels_transfer"][f"{target}_pred"]["method"])
        
        self.assertEqual(6, output_data.mod["rna"].uns["labels_transfer"][f"{targets[0]}_pred"]["max_depth"])
        self.assertEqual(6, output_data.mod["rna"].uns["labels_transfer"][f"{targets[1]}_pred"]["max_depth"])
        self.assertEqual(4, output_data.mod["rna"].uns["labels_transfer"][f"{targets[2]}_pred"]["max_depth"])

        # Remove output file to prevent errors in the next tests
        Path("output.h5mu").unlink()

if __name__ == '__main__':
    unittest.main()
