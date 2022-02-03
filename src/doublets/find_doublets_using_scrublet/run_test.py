import unittest
import os
from os import path
import subprocess
import scanpy as sc
import pandas
import numpy as np

class TestScrublet(unittest.TestCase):
    def test_simple_dedoubling(self):
        out = subprocess.check_output([
              "./find_doublets_using_scrublet", "--input", "pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5ad",  "--output=output-py1.h5ad"]
            ).decode("utf-8")
       
        self.assertTrue(path.exists("output-py1.h5ad"))
        
        data = sc.read_h5ad("output-py1.h5ad")
        self.assertTrue("ScrubletPredictedDoublets" in data.obs.columns)
        self.assertTrue("ScrubletDoubletScore" in data.obs.columns)

    def test_simple_dedoubling_csv(self):
        out = subprocess.check_output([
              "./find_doublets_using_scrublet", "--input", "pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5ad",  "--output=output-py2.csv", "--outputFormat", "csv"]
            ).decode("utf-8")
       
        self.assertTrue(path.exists("output-py2.csv"))
        
        data = pandas.read_csv("output-py2.csv")

        self.assertTrue("ScrubletPredictedDoublets" in data.columns)
        self.assertTrue("ScrubletDoubletScore" in data.columns)

    def test_simple_dedoubling_with_column_renaming(self):
        out = subprocess.check_output([
              "./find_doublets_using_scrublet", "--input", "pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5ad",  "--output=output-py2.csv", "--outputFormat", "csv", "--colNamePredictedDoublets", "scrublet-doublets"]
            ).decode("utf-8")
       
        self.assertTrue(path.exists("output-py2.csv"))
        
        data = pandas.read_csv("output-py2.csv")

        self.assertTrue("scrublet-doublets" in data.columns)
        self.assertTrue("ScrubletDoubletScore" in data.columns)
        


unittest.main()
