import subprocess
import mudata as mu
import logging
from sys import stdout
from unittest import TestCase, main
from pathlib import Path

## VIASH START
meta = {
    'functionality_name': './target/docker/filter/subset_h5mu/subset_h5mu',
    'resources_dir': 'resources_test/pbmc_1k_protein_v3/'
}
## VIASH END

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

resources_dir, functionality_name = meta["resources_dir"], meta["functionality_name"]
input_path = f"{meta['resources_dir']}/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"

class TestFilterWithCounts(TestCase):
    def setUp(self) -> None:
        mu_in = mu.read_h5mu(input_path)
        self.orig_obs = mu_in.mod['rna'].n_obs
        self.orig_vars = mu_in.mod['rna'].n_vars
        self.orig_prot_obs = mu_in.mod['prot'].n_obs
        self.orig_prot_vars = mu_in.mod['prot'].n_vars
        return super().setUp()

    def tearDown(self) -> None:
        return super().tearDown()

    def _run_and_check_output(self, args_as_list):
        try:
            subprocess_args = [meta['executable']] + args_as_list
            logger.info(" ".join(subprocess_args))
            subprocess.check_output(subprocess_args)
        except subprocess.CalledProcessError as e:
            logger.info(e.stdout.decode("utf-8"))
            raise e

    def test_filter_nothing(self):
        self._run_and_check_output([
            "--input", input_path,
            "--output", "output-1.h5mu",
            "--number_of_observations", "100",
            "--output_compression", "gzip"
            ])
        self.assertTrue(Path("output-1.h5mu").is_file(), msg="Output file not found")
        mu_out = mu.read_h5mu("output-1.h5mu")
        new_obs = mu_out.mod['rna'].n_obs
        new_vars = mu_out.mod['rna'].n_vars
        self.assertEqual(new_obs, 100, msg="No RNA obs should have been filtered")
        self.assertEqual(new_vars, self.orig_vars, msg="No RNA vars should have been filtered")
        self.assertEqual(mu_out.mod['prot'].n_obs, self.orig_prot_obs, msg="No prot obs should have been filtered")
        self.assertEqual(mu_out.mod['prot'].n_vars, self.orig_prot_vars, msg="No prot vars should have been filtered")
        self.assertListEqual(list(mu_out.mod['rna'].var['feature_types'].cat.categories), ["Gene Expression"],
                             msg="Feature types of RNA modality should be Gene Expression")
        self.assertListEqual(list(mu_out.mod['prot'].var['feature_types'].cat.categories), ["Antibody Capture"],
                             msg="Feature types of prot modality should be Antibody Capture")


if __name__ == "__main__":
    main()
