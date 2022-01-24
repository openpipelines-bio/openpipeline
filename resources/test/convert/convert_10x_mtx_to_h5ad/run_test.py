import unittest
import os
from os import path
import subprocess
import scanpy as sc
import pandas
import numpy as np

class TestMtxToH5AD(unittest.TestCase):
    def test_simple_conversion(self):
        out = subprocess.check_output([
              "./mtxtoh5ad", "--input", ".",  "--output=output-py1.h5ad"]
            ).decode("utf-8")
       
        self.assertTrue(path.exists("output-py1.h5ad"))
        
        data = sc.read_h5ad("output-py1.h5ad")
        self.assertEqual(len(data), 996)

unittest.main()