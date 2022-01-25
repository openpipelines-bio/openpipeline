import unittest
import os
from os import path
import subprocess
import scanpy as sc
import pandas

class TestFilter(unittest.TestCase):
    def test_check_gene_expression_filtering(self):
        initialData = sc.read_h5ad("pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5ad")
        self.assertEqual(len(initialData), 713)

        out = subprocess.check_output([
              "./filter_on_umi_genes_and_mito", 
              "--input", "pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5ad", "--output=output-1.h5ad", 
              ]
            ).decode("utf-8")

        self.assertTrue(path.exists("output-1.h5ad"))
        
        data = sc.read_h5ad("output-1.h5ad")

        self.assertCountEqual(data.var["feature_types"].cat.categories, ["Gene Expression"])
        self.assertEqual(len(data), 567)


    def test_check_gene_expression_filtering_2(self):
        initialData = sc.read_h5ad("pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5ad")
        self.assertEqual(len(initialData), 713)

        out = subprocess.check_output([
              "./filter_on_umi_genes_and_mito", 
              "--input", "pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5ad", "--output=output-1.h5ad", 
              "--minUMICount=2000"
              ]
            ).decode("utf-8")

        self.assertTrue(path.exists("output-1.h5ad"))
        
        data = sc.read_h5ad("output-1.h5ad")

        self.assertCountEqual(data.var["feature_types"].cat.categories, ["Gene Expression"])
        self.assertEqual(len(data), 548)



unittest.main()


    
