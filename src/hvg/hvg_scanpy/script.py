### VIASH START
par = {
    "input": "resources/test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "output": "output.h5mu",
    "excluded_genes": "",
    "flavor": "seurat",
}
### VIASH END

import scanpy as sc
import muon as mu

mdata = mu.read_h5mu(par["input"])

sc.pp.log1p(mdata.mod["rna"])
sc.pp.highly_variable_genes(mdata.mod["rna"], flavor=par["flavor"])

if len(par["excluded_genes"]) > 0:
    excluded_genes = list(map(str.strip, par["excluded_genes"].split(",")))
    print("Excluding genes: " + str(excluded_genes))

    mdata.var["highly_variable"] = (
        ~mdata.var["highly_variable"].index.isin(excluded_genes)
        & mdata.var["highly_variable"].values
    )

mdata.write_h5mu(filename=par["output"])
