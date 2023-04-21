from pathlib import Path
from unittest import main, TestCase
import subprocess
import mudata as mu
import logging
from sys import stdout
from tempfile import NamedTemporaryFile 
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
import re

## VIASH START
meta = {
    'functionality_name': 'foo',
    'resources_dir': 'resources_test/',
    'executable': '/home/di/code/openpipeline/target/docker/filter/filter_with_scrublet/filter_with_scrublet'
}
## VIASH END

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

resources_dir, functionality_name = meta["resources_dir"], meta["functionality_name"]
input_path = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"


class TestFilterWithScrublet(TestCase):
    def setUp(self) -> None:
        self.tempfile = NamedTemporaryFile(suffix=".h5mu")
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
            subprocess.check_output(subprocess_args, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            logger.info(e.stdout.decode("utf-8"))
            raise e

    def test_filter_a_little_bit(self):
        self._run_and_check_output([
                "--input", input_path,
                "--output", "output-1.h5mu",
                "--min_counts", "3",
                "--output_compression", "gzip"
        ])
        self.assertTrue(Path("output-1.h5mu").is_file(), msg="Output file not found")
        mu_out = mu.read_h5mu("output-1.h5mu")
        self.assertIn("filter_with_scrublet", mu_out.mod["rna"].obs)
        new_obs = mu_out.mod['rna'].n_obs
        new_vars = mu_out.mod['rna'].n_vars
        self.assertEqual(new_obs, self.orig_obs, msg="No RNA obs should have been filtered")
        self.assertEqual(new_vars, self.orig_vars, msg="No RNA vars should have been filtered")
        self.assertEqual(mu_out.mod['prot'].n_obs, self.orig_prot_obs, msg="No prot obs should have been filtered")
        self.assertEqual(mu_out.mod['prot'].n_vars, self.orig_prot_vars, msg="No prot vars should have been filtered")
        self.assertListEqual(list(mu_out.mod['rna'].var['feature_types'].cat.categories), ["Gene Expression"],
                             msg="Feature types of RNA modality should be Gene Expression")
        self.assertListEqual(list(mu_out.mod['prot'].var['feature_types'].cat.categories), ["Antibody Capture"],
                             msg="Feature types of prot modality should be Antibody Capture")

    def test_filtering_a_lot(self):
        self._run_and_check_output([
            "--input", input_path,
            "--output", "output-2.h5mu",
            "--modality", "rna",
            "--min_counts", "10",
            "--num_pca_components", "10",
            "--do_subset"
            ])
        self.assertTrue(Path("output-2.h5mu").is_file(), msg="Output file not found")
        mu_out = mu.read_h5mu("output-2.h5mu")
        new_obs = mu_out.mod['rna'].n_obs
        new_vars = mu_out.mod['rna'].n_vars
        self.assertLess(new_obs, self.orig_obs, msg="Some cells should have been filtered")
        self.assertEqual(new_vars, self.orig_vars, msg="No genes should have been filtered")
        self.assertEqual(mu_out.mod['prot'].n_obs, self.orig_obs, msg="No prot obs should have been filtered")
        self.assertEqual(mu_out.mod['prot'].n_vars, self.orig_prot_vars, msg="No prot vars should have been filtered")
        self.assertListEqual(list(mu_out.mod['rna'].var['feature_types'].cat.categories), ["Gene Expression"],
                             msg="Feature types of RNA modality should be Gene Expression")
        self.assertListEqual(list(mu_out.mod['prot'].var['feature_types'].cat.categories), ["Antibody Capture"],
                             msg="Feature types of prot modality should be Antibody Capture" )

    def test_doublet_automatic_threshold_detection_fails_recovery(self):
        """
        Test if the component can recover from scrublet not automatically able to set the doublet score threshold
        and it is not set.
        """
        np.random.seed(4) # Make test reproducable
        mudata = mu.read_h5mu(input_path)
        # Simulate a failed scrublet run by passing very little cells
        mudata = mudata[152].copy()
        nobs = 100
        x_data = np.repeat(mudata.mod['rna'].X.todense(), nobs, axis=0)
        replace_rate = 0.000001
        # Random perturbations because otherwise the detection fails in other ways (PCA cannot be run)
        mask = np.random.choice([0, 1], size=x_data.shape, p=((1 - replace_rate), replace_rate)).astype("bool")
        r = np.random.rand(*x_data.shape)*np.max(x_data)
        x_data[mask] = r[mask]
        obs_name = mudata.mod['rna'].obs.index.to_list()[0]
        obs_data = pd.DataFrame([], index=[f"{obs_name}_{i}" for i in range(nobs)])
        mod = ad.AnnData(X=csr_matrix(x_data), obs=obs_data,var=mudata.mod['rna'].var)
        new_mudata = mu.MuData({"rna": mod})
        new_mudata.update()
        with NamedTemporaryFile(suffix="h5mu") as temp_file:
            new_mudata.write(temp_file.name)
            try:
                self._run_and_check_output([
                    "--input", temp_file.name,
                    "--output", "output-4.h5mu",
                    "--output_compression", "gzip",
                    "--num_pca_components", "1",
                    "--min_gene_variablity_percent", "0"
                ])
            except subprocess.CalledProcessError as e:
                assert re.search(r"RuntimeError: Scrublet could not automatically detect the doublet score threshold\. "
                                 r"--allow_automatic_threshold_detection_fail can be used to ignore this failure and "
                                r"set the corresponding output columns to NA\.",
                        e.stdout.decode('utf-8'))

            self._run_and_check_output([
                "--input", temp_file.name,
                "--output", "output-5.h5mu",
                "--output_compression", "gzip",
                "--num_pca_components", "1",
                "--min_gene_variablity_percent", "0",
                "--allow_automatic_threshold_detection_fail"
            ])
            self.assertTrue(Path("output-5.h5mu").is_file(), msg="Output file not found")
            mu_out = mu.read_h5mu("output-5.h5mu")
            assert not mu_out.mod['rna'].obs['filter_with_scrublet'].any()



if __name__ == "__main__":
    main()
