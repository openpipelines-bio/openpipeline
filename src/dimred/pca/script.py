import scanpy as sc
import muon as mu

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "output": "output.h5mu",
    "nPCs": 25,
}
## VIASH END

print("Reading", par["input"])
mdata = mu.read_h5mu(par["input"])

print("Computing PCA coordinates")
sc.tl.pca(mdata.mod["rna"], n_comps=par["nPCs"])
mdata.uns["pcaParameters"] = {"PCA: nPCs": par["nPCs"]}

print("Writing", par["output"])
mdata.write_h5mu(filename=par["output"])
