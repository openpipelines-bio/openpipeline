import mudata as mu

##VIASH START
par = {
    "input": "foo.final.h5mu",
}

meta = {"resources_dir": "resources_test/pbmc_1k_protein_v3"}

##VIASH END

print("Loading data", flush=True)
data = mu.read_h5mu(par["input"])

assert "X_umap" in data.mod["rna"].obsm, "X_umap not found in .obsm"
assert data.mod["rna"].obsm["X_umap"].shape[1] == 2, (
    f"X_umap has wrong shape expected 2 n_comp but got {data.mod['rna'].obsm['X_umap'].shape[1]}"
)
assert "pca_variance" in data.mod["rna"].uns
assert "pca_loadings" in data.mod["rna"].varm


print("Test successful!", flush=True)
