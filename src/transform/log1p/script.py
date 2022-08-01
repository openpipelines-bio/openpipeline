import scanpy as sc
import muon as mu

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "output": "output.h5mu",
    "base": None,
    "modality": ["rna"],
}
meta = {"functionality_name": "lognorm"}
## VIASH END

print("Reading input mudata")
mdata = mu.read_h5mu(par["input"])
mdata.var_names_make_unique()


for mod in par["modality"]:
    print(f"Performing log transformation on modality {mod}")
    data = mdata.mod[mod]
    new_layer = sc.pp.log1p(data,
                            base=par["base"],
                            copy=True if par['output_layer'] else False)
    if new_layer:
        data.layers[par['output_layer']] = new_layer.X

print("Writing to file")
mdata.write(filename=par["output"])
