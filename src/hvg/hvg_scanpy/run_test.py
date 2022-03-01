import unittest
import os
from os import path
import subprocess
import scanpy as sc
import pandas
import numpy as np


class TestHVGSelection(unittest.TestCase):
    def test_simple_hvg(self):
        inputData = sc.read_h5ad("pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5ad")
        sc.pp.log1p(inputData)
        inputData.write_h5ad("lognormed.h5ad")

        self.assertTrue("highly_variable" not in inputData.var.columns)

        out = subprocess.check_output(
            ["./hvg_scanpy", "--input", "lognormed.h5ad", "--output=output-py1.h5ad"]
        ).decode("utf-8")

        self.assertTrue(path.exists("output-py1.h5ad"))
        data = sc.read_h5ad("output-py1.h5ad")

        self.assertTrue("highly_variable" in data.var.columns)

    def test_simple_hvg_with_exclusion(self):
        inputData = sc.read_h5ad("pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5ad")
        sc.pp.log1p(inputData)
        inputData.write_h5ad("lognormed.h5ad")

        self.assertTrue("highly_variable" not in inputData.var.columns)

        out = subprocess.check_output(
            [
                "./hvg_scanpy",
                "--input",
                "lognormed.h5ad",
                "--output=output-py1.h5ad",
                "--excluded_genes",
                "GNAI2,RHOA",
            ]
        ).decode("utf-8")

        self.assertTrue(path.exists("output-py1.h5ad"))
        data = sc.read_h5ad("output-py1.h5ad")

        self.assertTrue("highly_variable" in data.var.columns)
        self.assertFalse(data[:, data.var_names == "GNAI2"].var["highly_variable"][0])
        self.assertFalse(data[:, data.var_names == "RHOA"].var["highly_variable"][0])

        self.assertFalse(data[:, data.var_names == "JUN"].var["highly_variable"][0])


unittest.main()
