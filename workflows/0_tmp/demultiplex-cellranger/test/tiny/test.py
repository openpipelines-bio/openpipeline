import unittest
import os
from os import path
import subprocess
import scanpy as sc
import pandas
import numpy as np
import sys

class TestBCLGeneration(unittest.TestCase):
    def test_output_files_presence(self):
        self.assertTrue(path.isdir(self.INPUT_DIRECTORY), "The output directory was not found.")

        self.assertTrue(path.isdir(self.INPUT_DIRECTORY + "/fastq"), "No output fastq could be found in the output.")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        TestBCLGeneration.INPUT_DIRECTORY = sys.argv.pop()
        unittest.main()
