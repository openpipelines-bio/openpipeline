

### VIASH START

par = {
	"input": "./test/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5",
	"output": "./test/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5ad",
	"compression": "gzip",
	"conversions-obsm": "{'counts_antibody': 'adt'}"
}
### VIASH END

import muon as mu
import anndata
import json

data =  anndata.read_h5ad(par["input"])       

try: 
	data.var_names_make_unique()
except:
	pass

muon = mu.MuData({"rna": data})

for key, value in json.loads(par["conversions_obsm"]).items():
	if key in data.obsm:
		muon.mod[value] = anndata.AnnData(data.obsm[key])
		del muon["rna"].obsm[key]


muon.write_h5mu(par["output"], compression = par["compression"])
