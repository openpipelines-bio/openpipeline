import unittest
import os
from os import path
import subprocess
import scanpy as sc
import muon as mu
import pandas
import numpy as np

class TestH5toH5AD(unittest.TestCase):
    def test_simple_conversion(self):
        preData = sc.read_h5ad("pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5.h5ad")

        preData.obsm["counts_antibody"] = preData.to_df().iloc[:,0:10]
        preData.obsm["counts_crispr"] = preData.to_df().iloc[:,0:10]
        preData.obsm["counts_custom"] = preData.to_df().iloc[:,0:10]

        preData.write_h5ad("input-1.h5ad")

        out = subprocess.check_output([
              "./convert_h5ad_to_h5mu", 
              "--input=input-1.h5ad",  
              "--output=output-py1-rna.h5mu"]
            ).decode("utf-8")
       
        self.assertTrue(path.exists("output-py1-rna.h5mu"))
        
        data = mu.read_h5mu("output-py1-rna.h5mu")

        self.assertTrue("rna" in data.mod)
        self.assertTrue("prot" in data.mod)
        self.assertTrue("custom" in data.mod)

        self.assertFalse("counts_antibody" in data.mod["rna"].obsm)
        self.assertFalse("counts_custom" in data.mod["rna"].obsm)
        self.assertTrue("counts_crispr" in data.mod["rna"].obsm)


unittest.main()
