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

assert "rna" in input_mudata.mod.keys()
assert input_mudata.n_obs == 713
assert input_mudata.mod["rna"].var["feature_types"].unique() == [
    "Gene Expression"
], "Output X should only contain Gene Expression vars."
    