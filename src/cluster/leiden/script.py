import muon as mu
import scanpy as sc

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.norm.hvg.pca.nn.umap.h5mu",
    "output": "output.h5mu",
    "output_format": "h5mu",
    "cluster_column_name": "leiden.res.0.25",
    "resolution": float("0.25"),
}
## VIASH END

print("Reading", par["input"])
mdata = mu.read_h5mu(par["input"])

print("Cluster cells using the Leiden algorithm")
sc.tl.leiden(
    mdata.mod["rna"],
    resolution=float(par["resolution"]),
    key_added=par["cluster_column_name"],
)

print("Writing", par["output"])
mdata.write_h5mu(filename=par["output"])
