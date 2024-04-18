from mudata import read_h5mu

##VIASH START
par = {
    "input": "input.h5mu"
}

meta = {
    "resources_dir": "resources_test"
}
##VIASH END

input_mudata = read_h5mu(par["input"])
expected_colnames = ['gene_symbol', 'feature_types', 'genome']

assert list(input_mudata.mod.keys()) == ["rna"], "Input should contain rna modality."
assert list(input_mudata.var.columns) == expected_colnames, f"Input var columns should be: {expected_colnames}."
assert list(input_mudata.mod["rna"].var.columns) == expected_colnames, f"Input mod['rna'] var columns should be: {expected_colnames}."
