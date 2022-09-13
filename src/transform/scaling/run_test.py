import subprocess
import sys
import unittest
from pathlib import Path
from unittest import TestCase

import numpy as np
from mudata import read_h5mu

## VIASH START
meta = {
    'functionality_name': './target/native/transform/scale/scale',
    'resources_dir': './resources_test/pbmc_1k_protein_v3/'
}
## VIASH END

resources_dir, functionality_name = meta["resources_dir"], meta["functionality_name"]
input_file = f"{resources_dir}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"

class TestScaling(TestCase):
    def _run_and_check_output(self, args_as_list):
        try:
            subprocess.check_output([meta['executable']] + args_as_list)
        except subprocess.CalledProcessError as e:
            print(e.stdout.decode("utf-8"))
            raise e

    def test_scaling(self):
        """
        Output data must be centered around mean 0 and it has unit variance.
        """
        self._run_and_check_output([
            "--input", input_file,
            "--output", "scaled.h5mu"])

        self.assertTrue(Path("scaled.h5mu").is_file())
        output_data = read_h5mu("scaled.h5mu")
        output_x = output_data.mod['rna'].X
        mean = np.mean(output_x, axis=0, dtype=np.float64)
        variance = np.multiply(output_x, output_x).mean(axis=0, dtype=np.float64) - mean**2
        variance[variance == 0] = 1
        self.assertTrue(np.all(np.isclose(mean, 0, rtol=1e-07, atol=1e-07)))
        self.assertTrue(np.all(np.isclose(variance, 1, rtol=1e-03, atol=1e-03)))

    def test_scaling_noncenter(self):
        """
        Check if centering can be disabled.
        """
        self._run_and_check_output([
            "--input", input_file,
            "--output", "scaled.h5mu",
            "--zero_center", "false"])
        self.assertTrue(Path("scaled.h5mu").is_file())
        output_data = read_h5mu("scaled.h5mu")
        output_x = output_data.mod['rna'].X
        mean = np.mean(output_x, axis=0, dtype=np.float64)
        self.assertFalse(np.all(np.isclose(mean, 0, rtol=1e-07, atol=1e-07)))

    def test_scaling_maxvalue(self):
        """
        Check if output data is clipped when using --max_value
        """
        self._run_and_check_output([
            "--input", input_file,
            "--output", "scaled.h5mu",
            "--max_value", "0.5"])
        self.assertTrue(Path("scaled.h5mu").is_file())
        output_data = read_h5mu("scaled.h5mu")
        output_x = output_data.mod['rna'].X
        self.assertTrue(np.all(output_x <= 0.5))

    def test_scaling_modality(self):
        """
        Check if 'rna' modality remain untouched when using '--modality prot' argument.
        """
        self._run_and_check_output([
            "--input", input_file,
            "--output", "scaled.h5mu",
            "--modality", "prot"])
        self.assertTrue(Path("scaled.h5mu").is_file())
        input_data =  read_h5mu(input_file)
        output_data = read_h5mu("scaled.h5mu")
        output_rna = output_data.mod['rna'].X
        self.assertTrue(np.allclose(input_data.mod['rna'].X.todense(),
                                    output_rna.todense(), equal_nan=True))

        output_prot =  output_data.mod['prot'].X
        mean = np.mean(output_prot, axis=0, dtype=np.float64)
        variance = np.multiply(output_prot, output_prot).mean(axis=0, dtype=np.float64) - mean**2
        variance[variance == 0] = 1
        self.assertTrue(np.all(np.isclose(mean, 0, rtol=1e-07, atol=1e-07)))
        self.assertTrue(np.all(np.isclose(variance, 1, rtol=1e-03, atol=1e-03)))

if __name__ == "__main__":
    unittest.main()
