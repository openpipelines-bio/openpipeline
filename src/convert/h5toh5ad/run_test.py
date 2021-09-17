import unittest
import os
from os import path
import subprocess
import scanpy as sc
import pandas
import numpy as np

class TestH5toH5AD(unittest.TestCase):
    def test_simple_conversion(self):
        out = subprocess.check_output([
              "./h5toh5ad", "--input", "pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5",  "--output=output-py1.h5ad"]
            ).decode("utf-8")
       
        self.assertTrue(path.exists("output-py1.h5ad"))
        
        data = sc.read_h5ad("output-py1.h5ad")

        self.assertEqual(set(["Gene Expression"]), set(data.var["feature_types"].unique()))

        self.assertTrue("counts_antibody" in data.obsm)
        self.assertTrue("CD3_TotalSeqB" in data.obsm["counts_antibody"].columns)

unittest.main()
