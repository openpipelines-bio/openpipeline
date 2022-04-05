import muon as mu
import scanpy as sc

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.neighbors.h5mu",
    "output": "output.h5mu",
    "modality": ["rna"],
    "output_format": "h5mu",
    "obs_name": "leiden",
    "resolution": 0.25,
}
## VIASH END

print("Reading", par["input"])
mdata = mu.read_h5mu(par["input"])

for mod in par['modality']:
    print(f"Processing modality '{mod}'")
    data = mdata.mod[mod]
    sc.tl.leiden(
        data,
        resolution=par["resolution"],
        key_added=par["obs_name"],
    )

print("Writing", par["output"])
mdata.write_h5mu(filename=par["output"])
