import scanpy as sc

### VIASH START
par = {
    "input": "resources/test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.norm.hvg.pca.nn.umap.h5ad",
    "output": "output.h5ad",
    "output_format": "h5ad",
    "cluster_column_name": "leiden.res.0.25",
    "resolution": float("0.25"),
}
### VIASH END

print("Reading", par["input"])
data = sc.read_h5ad(par["input"])

print("Clustering cells using Leiden algorithm")
sc.tl.leiden(
    data, resolution=float(par["resolution"]), key_added=par["cluster_column_name"]
)

print("Writing file as", par["output_format"], "format", "to", par["output"])
if par["output_format"] == "h5ad":
    data.write_h5ad(par["output"], compression="lzf")
elif par["output_format"] == "csv":
    data.obs[par["cluster_column_name"]].to_csv(par["output"])
else:
    raise ValueError("An unrecognized output format was specified.")
