from mudata import read_h5mu
import bbknn

### VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "modality": "rna",
    "obs_batch": "sample_id",
    "obsm_input": "X_pca",
    "n_neighbors_within_batch": 3,
    "n_trim": None,
    "n_pcs": 50,
    "output": "output.h5mu",
    "output_compression": "gzip",
    "obsp_connectivities": "my_connectivities",
    "obsp_distances": "my_distances",
    "uns_output": "my_neighbors",
}
### VIASH END

mudata = read_h5mu(par["input"])
adata = mudata.mod[par["modality"]]

# copy data
tmp_adata = adata.copy()
bbknn.bbknn(
    tmp_adata,
    use_rep=par["obsm_input"],
    batch_key=par["obs_batch"],
    neighbors_within_batch=par["n_neighbors_within_batch"],
    n_pcs=par["n_pcs"],
    trim=par["n_trim"],
)

# store output
adata.obsp[par["obsp_connectivities"]] = tmp_adata.obsp["connectivities"]
adata.obsp[par["obsp_distances"]] = tmp_adata.obsp["distances"]
adata.uns[par["uns_output"]] = tmp_adata.uns["neighbors"]
adata.uns[par["uns_output"]]["distances_key"] = par["obsp_distances"]
adata.uns[par["uns_output"]]["connectivities_key"] = par["obsp_connectivities"]

# write to file
mudata.write_h5mu(par["output"], compression=par["output_compression"])
