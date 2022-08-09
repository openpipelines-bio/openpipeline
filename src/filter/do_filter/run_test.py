from pathlib import Path
import subprocess
import muon
import numpy as np
from unittest import TestCase, main
from tempfile import NamedTemporaryFile
## VIASH START
meta = {
    'functionality_name': './target/native/filter/do_filter/do_filter',
    'resources_dir': 'resources_test/'
}
## VIASH END
resources_dir, functionality_name = meta["resources_dir"], meta["functionality_name"]
input_path = f"{resources_dir}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"


class TestDoFilter(TestCase):
    def setUp(self) -> None:
        self._tempfile = NamedTemporaryFile(suffix=".h5mu")
        mu_in = muon.read_h5mu(input_path)
        self.orig_obs = mu_in.mod['rna'].n_obs
        self.orig_vars = mu_in.mod['rna'].n_vars

        ad_rna = mu_in.mod['rna']
        print(f"  input: {ad_rna}")
        ad_rna.obs["filter_none"] = np.repeat(True, ad_rna.n_obs)
        ad_rna.obs["filter_with_random"] = np.random.choice(a=[False, True], size=ad_rna.n_obs)
        ad_rna.var["filter_with_random"] = np.random.choice(a=[False, True], size=ad_rna.n_vars) 
        mu_in.write_h5mu(self._tempfile.name)       
        return super().setUp()
    
    def tearDown(self) -> None:
        self._tempfile.close()
        return super().tearDown()

    def _run_and_check_output(self, args_as_list):
        try:
            subprocess_args = [f"./{functionality_name}"] + args_as_list
            print(" ".join(subprocess_args))
            subprocess.check_output(subprocess_args)
        except subprocess.CalledProcessError as e:
            print(e.stdout.decode("utf-8"))
            raise e
    
    def test_filtering_a_little_bit(self):
        self._run_and_check_output([
            "--input", self._tempfile.name, 
            "--output", "output-1.h5mu",
            "--obs_filter", "filter_none:filter_with_random",
            "--var_filter", "filter_with_random"]
        )
        self.assertTrue(Path("output-1.h5mu").is_file(), msg="Output file not found")
        mu_out = muon.read_h5mu("output-1.h5mu")
        print(f"  output1: {mu_out.mod['rna']}")
        new_obs = mu_out.mod['rna'].n_obs
        new_vars = mu_out.mod['rna'].n_vars
        self.assertTrue(new_obs < self.orig_obs, msg="Some RNA obs should have been filtered")
        self.assertTrue(new_vars < self.orig_vars, msg="Some RNA vars should have been filtered")

    def test_filter_nothing(self):
        self._run_and_check_output([
            "--input", self._tempfile.name, 
            "--output", "output-2.h5mu",
            "--obs_filter", "filter_none"])
        self.assertTrue(Path("output-2.h5mu").is_file(), msg="Output file not found")
        mu_out = muon.read_h5mu("output-2.h5mu")
        print(f"  output2: {mu_out.mod['rna']}")
        new_obs = mu_out.mod['rna'].n_obs
        new_vars = mu_out.mod['rna'].n_vars
        self.assertEqual(new_obs, self.orig_obs, msg="No RNA obs should have been filtered")
        self.assertEqual(new_vars, self.orig_vars, msg="No RNA vars should have been filtered")

if __name__ == "__main__":
    main()