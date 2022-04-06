import muon as mu
import scanpy as sc

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.tx_processing.h5mu",
    "output": "output.h5mu",
    "metric": 'cosine',
    "num_neighbors": 15,
    "modality": ["rna"],
    "obsp_name_prefix": "foo"
}
meta = {"functionality_name": "lognorm"}
## VIASH END

print("Reading input mudata")
mdata = mu.read_h5mu(par["input"])
mdata.var_names_make_unique()

for mod in par["modality"]:
    print(f"Computing a neighborhood graph on modality {mod}")
    sc.pp.neighbors(
        mdata.mod[mod],
        n_neighbors=par["num_neighbors"], 
        metric=par["metric"],
        key_added=par["obsp_name_prefix"]
    )

print("Writing to file")
mdata.write(filename=par["output"])
