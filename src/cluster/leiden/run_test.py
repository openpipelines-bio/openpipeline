import unittest
import os
from os import path
import subprocess
import scanpy as sc
import pandas
import numpy as np

class TestLeiden(unittest.TestCase):
    def test_simple_clustering(self):
        out = subprocess.check_output([
              "./leiden", "--input", "pbmc_1k_protein_v3_filtered_feature_bc_matrix.norm.hvg.pca.nn.umap.h5ad",  "--output=output-py1.h5ad"]
            ).decode("utf-8")
       
        self.assertTrue(path.exists("output-py1.h5ad"))
        
        data = sc.read_h5ad("output-py1.h5ad")
        self.assertTrue("leiden.res.0.25" in data.obs.columns)
        
    def test_simple_clustering_with_resolution_1(self):
         out = subprocess.check_output([
              "./leiden", "--input", "pbmc_1k_protein_v3_filtered_feature_bc_matrix.norm.hvg.pca.nn.umap.h5ad",  "--output=output-py2.h5ad", "--resolution=1", "--clusterColumnName=leiden.custom.resolution"]
            ).decode("utf-8")

         self.assertTrue(path.exists("output-py2.h5ad"))

         data = sc.read_h5ad("output-py2.h5ad")
         self.assertTrue("leiden.custom.resolution" in data.obs.columns)

    def test_csv_output(self):
         out = subprocess.check_output([
              "./leiden", "--input", "pbmc_1k_protein_v3_filtered_feature_bc_matrix.norm.hvg.pca.nn.umap.h5ad",  "--output=output-py3.csv", "--outputFormat=csv", "--resolution=1", "--clusterColumnName=leiden.custom.resolution"]
            ).decode("utf-8")

         self.assertTrue(path.exists("output-py3.csv"))
         data = pandas.read_csv("output-py3.csv")
         self.assertTrue("leiden.custom.resolution" in data.columns) 


unittest.main()
