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

assert list(input_mudata.mod.keys()) == ['rna', 'prot', 'vdj_t']
assert list(input_mudata.uns.keys()) == ['metrics_cellranger']
expected_metrics = ['Category', 'Library Type', 'Grouped By', 'Group Name', 'Metric Name', 'Metric Value']
assert input_mudata.uns['metrics_cellranger'].columns.to_list() == expected_metrics
