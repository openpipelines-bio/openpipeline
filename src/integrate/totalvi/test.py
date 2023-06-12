from pathlib import Path
import subprocess
import unittest

import mudata

## VIASH START
meta = {
    "executable": "./target/docker/integrate/totalvi/totalvi",
    "resources_dir": "./resources_test/pbmc_1k_protein_v3/"
}
## VIASH END

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3_mms.h5mu"

class TesttotalVI(unittest.TestCase):
    def _run_and_check_output(self, args_as_list):
        try:
            subprocess.check_output([meta["executable"]] + args_as_list)
        except subprocess.CalledProcessError as e:
            print(e.stdout.decode("utf-8"))
            raise e

    def test_totalvi(self):
        """Map data containing proteins on itself"""
        self._run_and_check_output([
            "--input", input_file,
            "--reference", input_file,
            "--query_proteins_modality", "prot",
            "--reference_proteins_modality", "prot",
            "--var_input", "filter_with_hvg",
            "--reference_model_path", "totalvi_reference_model/",
            "--query_model_path", "totalvi_query_model/",
            "--max_epochs", "1",
            "--max_query_epochs", "1",
            "--output", "output.h5mu"])

        self.assertTrue(Path("output.h5mu").is_file())
        output_data = mudata.read_h5mu("output.h5mu")
        self.assertIn("X_integrated_totalvi", output_data.mod["rna"].obsm)
        self.assertIn("_scvi_batch", output_data.mod["rna"].obs.columns.tolist())
        self.assertIn("_scvi_labels",  output_data.mod["rna"].obs.columns.tolist())
        self.assertIn("X_totalvi_normalized_rna", output_data.mod["rna"].obsm)
        self.assertIn("X_totalvi_normalized_protein", output_data.mod["rna"].obsm)


if __name__ == "__main__":
    unittest.main()