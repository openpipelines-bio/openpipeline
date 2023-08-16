import unittest
import subprocess
from tempfile import NamedTemporaryFile
import mudata
from pathlib import Path
import numpy as np

## VIASH START
meta = {
    'executable': './target/native/interpret/lianapy/',
    'resources_dir': './resources_test/pbmc_1k_protein_v3/'
}
## VIASH END

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3_mms.h5mu"

class TestLianaPy(unittest.TestCase):
    def _run_and_check_output(self, args_as_list):
        try:
            subprocess.check_output([meta['executable']] + args_as_list)
        except subprocess.CalledProcessError as e:
            print(e.stdout.decode("utf-8"))
            raise e

    def test_lianapy(self):
        input_data = mudata.read_h5mu(input_file)
        self._run_and_check_output([
            "--input", input_file,
            "--output_compression", "gzip",
            "--modality", "rna",
            # "--layer", "log_normalized",
            "--groupby", "bulk_labels",
            "--resource_name", "consensus",
            "--gene_symbol", "gene_symbol",
            "--expr_prop", "0.1",
            "--min_cells", "5",
            "--aggregate_method", "rra",
            "--return_all_lrs", "False",
            "--n_perms", "11",
            "--output", "output.h5mu"])
        self.assertTrue(Path("output.h5mu").is_file())
        output_data = mudata.read_h5mu("output.h5mu")
        np.testing.assert_array_equal(output_data.mod['rna'].X.data, input_data.mod['rna'].X.data)
        np.testing.assert_array_equal(input_data.mod['rna'].var.index, output_data.mod['rna'].var.index)
        assert "liana_res" in output_data.mod["rna"].uns
        all(elem in output_data.mod['rna'].obs['bulk_labels'].values for elem in output_data.mod['rna'].uns['liana_res']['source'].unique())
        all(elem in output_data.mod['rna'].obs['bulk_labels'].values for elem in output_data.mod['rna'].uns['liana_res']['target'].unique())

if __name__ == '__main__':
    unittest.main()