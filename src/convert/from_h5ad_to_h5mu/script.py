import muon as mu
import anndata
import json

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5ad",
    "modality": "rna",
    "output": "output.h5mu",
    "compression": "gzip",
    "conversions_obsm": '{"counts_antibody":"prot", "counts_custom": "custom"}',
}
## VIASH END

print("Reading", par["input"])
data = anndata.read_h5ad(par["input"])

try:
    data.var_names_make_unique()
except:
    pass

print("Converting ")
muon = mu.MuData({par["rna"]: data})

for source, target in json.loads(par["conversions_obsm"]).items():
    if source in data.obsm:
        muon.mod[target] = anndata.AnnData(data.obsm[source])
        del muon[par["rna"]].obsm[source]

print(f"Writing {par['output']}")
muon.write_h5mu(par["output"])
