### VIASH START
par = {
    "input": "resources/test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "output": "output.h5mu",
    "min_dist": 0.5,
    "alpha": 1.0,
    "gamma": 1.0,
    "random_seed": 0,
    "negative_sample_rate": 5,
    "init_pos": "spectral",
}
### VIASH END

import scanpy as sc
import muon as mu

print("Reading", par["input"])
mdata = mu.read_h5mu(par["input"])

# TODO Note by Eric: Not sure if this is called earlier in the pipeline already? UMAP won't work on the h5mu in the resources folder without this.
print("Computing a neighborhood graph of observations")
sc.pp.neighbors(mdata.mod["rna"], use_rep="X")

print(
    "Embedding the neighborhood graph using Uniform Manifold Approximation and Projection (UMAP)"
)

sc.tl.umap(
    mdata.mod["rna"],
    min_dist=par["min_dist"],
    alpha=par["alpha"],
    gamma=par["gamma"],
    random_state=par["random_seed"],
    negative_sample_rate=par["negative_sample_rate"],
    init_pos=par["init_pos"],
)

print("Writing", par["output"])
mdata.write_h5mu(filename=par["output"])
