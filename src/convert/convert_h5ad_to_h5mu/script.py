import muon as mu
import anndata
import json

### VIASH START
par = {
    "input": "resources/test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5ad",
    "output": "output.h5mu",
    "compression": "gzip",
    "conversions_obsm": '{"counts_antibody":"prot", "counts_custom": "custom"}',
}
### VIASH END

data = anndata.read_h5ad(par["input"])

try:
    data.var_names_make_unique()
except:
    pass

muon = mu.MuData({"rna": data})

for key, value in json.loads(par["conversions_obsm"]).items():
    if key in data.obsm:
        muon.mod[value] = anndata.AnnData(data.obsm[key])
        del muon["rna"].obsm[key]

muon.write_h5mu(par["output"], compression=par["compression"])
