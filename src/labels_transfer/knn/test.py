import unittest
import subprocess
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

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"
reference_file = f"{meta['resources_dir']}/annotation_test_data/TS_Blood_filtered.h5ad"

class TestKNN(unittest.TestCase):
    def _run_and_check_output(self, args_as_list):
        try:
            subprocess.check_output([meta['executable']] + args_as_list)
        except subprocess.CalledProcessError as e:
            print(e.stdout.decode("utf-8"))
            raise e

    def test_one_class(self):
        with NamedTemporaryFile("w", suffix=".h5mu") as tempfile_reference_file:
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

                self._run_and_check_output([
                    "--input", tempfile_input_file.name,
                    "--modality", "rna",
                    "--query_obsm_key", "X_integrated_scanvi",
                    "--reference", tempfile_reference_file.name,
                    "--reference_obsm_key", "X_integrated_scanvi",
                    "--output", "output.h5mu",
                    "-k", "5",
                    "--targets", "ann_level_1"])

        self.assertTrue(Path("output.h5mu").is_file())

        output_data = mudata.read_h5mu("output.h5mu")
        self.assertIn("ann_level_1_pred", output_data.mod["rna"].obs)
        self.assertIn("ann_level_1_uncertainty", output_data.mod["rna"].obs)
        self.assertIn("labels_transfer", output_data.mod["rna"].uns)
        self.assertIn("ann_level_1_pred", output_data.mod["rna"].uns["labels_transfer"])
        self.assertEqual("KNN_pynndescent", output_data.mod["rna"].uns["labels_transfer"]["ann_level_1_pred"]["method"])
        self.assertEqual(5, output_data.mod["rna"].uns["labels_transfer"]["ann_level_1_pred"]["n_neighbors"])
        self.assertNotIn("ann_level_2_pred", output_data.mod["rna"].obs)

    def test_multiple_classes(self):

        targets = ["ann_level_1", "ann_level_2", "ann_level_3", "ann_level_4", "ann_level_5", "ann_finest_level"]
        
        with NamedTemporaryFile("w", suffix=".h5mu") as tempfile_reference_file:
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
                    "--output", "output.h5mu",
                    "-k", "5",
                    "--targets", ",".join(targets)])

        self.assertTrue(Path("output.h5mu").is_file())

        output_data = mudata.read_h5mu("output.h5mu")

        for target in targets:
            self.assertIn(f"{target}_pred", output_data.mod["rna"].obs)
            self.assertIn(f"{target}_uncertainty", output_data.mod["rna"].obs)
            self.assertIn("labels_transfer", output_data.mod["rna"].uns)
            self.assertIn(f"{target}_pred", output_data.mod["rna"].uns["labels_transfer"])
            self.assertEqual("KNN_pynndescent", output_data.mod["rna"].uns["labels_transfer"][f"{target}_pred"]["method"])
            self.assertEqual(5, output_data.mod["rna"].uns["labels_transfer"][f"{target}_pred"]["n_neighbors"])

if __name__ == '__main__':
    unittest.main()
