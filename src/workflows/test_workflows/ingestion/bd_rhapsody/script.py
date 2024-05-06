from mudata import read_h5mu
import numpy as np

##VIASH START
par = {
    "input": "input.h5mu"
}

meta = {
    "resources_dir": "resources_test"
}
##VIASH END

input_mudata = read_h5mu(par["input"])
expected_var = ['gene_name', 'feature_types', 'reference_file']
expected_obs = ['run_id', 'library_id', 'sample_id']

assert list(input_mudata.mod.keys()) == ["rna"], "Input should contain rna modality."
assert list(input_mudata.var.columns) == expected_var, f"Input var columns should be: {expected_var}."
assert list(input_mudata.mod["rna"].var.columns) == expected_var, f"Input mod['rna'] var columns should be: {expected_var}."
assert list(input_mudata.mod["rna"].obs.columns) == expected_obs, f"Input obs columns should be: {expected_obs}."

assert np.array_equal(input_mudata.var["feature_types"].unique(), ["Gene Expression"]), "Output X should only contain Gene Expression vars."
