from mudata import read_h5mu, MuData


### VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "modality": ["rna"],
    "output": "foo.h5mu",
}
### VIASH END


input_mudata = read_h5mu(par["input"])
new_mods = {
    mod_name: mod
    for mod_name, mod in input_mudata.mod.items()
    if mod_name not in par["modality"]
}

new_mudata = MuData(new_mods)
new_mudata.write_h5mu(filename=par["output"], compression=par["output_compression"])
