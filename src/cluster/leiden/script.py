import scanpy as sc

### VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.norm.hvg.pca.nn.umap.h5mu",
    "output": "output.h5mu",
    "output_format": "h5mu",
    "cluster_column_name": "leiden.res.0.25",
    "resolution": float("0.25"),
}
### VIASH END

import scanpy as sc
import muon as mu

print("Reading", par["input"])
mdata = mu.read_h5mu(par["input"])

sc.tl.leiden(
    mdata.mod["rna"],
    resolution=float(par["resolution"]),
    key_added=par["cluster_column_name"],
)

mdata.write_h5mu(filename=par["output"])
