import json
import scanpy as sc
import scipy
import mudata as mu
import anndata
import logging
from sys import stdout

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.csv",
    "output": "output.h5mu",
    "delimiter": ",",
    "use_column_name": "true",
    "compression": "gzip",
}
## VIASH END
logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)


logger.info("Reading %s.", par["input"])
data = sc.read_csv(
    par["input"], delimiter=par["delimiter"], first_column_names=par["use_column_names"]
)

logger.info("Renaming keys.")
for adata in data.mod.values():
    adata.var.rename(columns={'gene_ids': 'gene_id', 'feature_types': 'feature_type'}, inplace=True)
data.var = data.var.drop(["feature_types", "gene_ids"], axis=1)
data.update()

logger.info("Converting.")
data.var_names_make_unique()
data.X = scipy.sparse.csr_matrix(data.X)
data.raw = data

mudata_obj = mu.MuData({"rna": data})

for key, value in json.loads(par["conversions_obsm"]).items():
    if key in data.obsm:
        mudata_obj.mod[value] = anndata.AnnData(data.obsm[key])
        del mudata_obj["rna"].obsm[key]

logger.info("Writing %s.", par["output"])
mudata_obj.write_h5mu(filename=par["output"])
logger.info("Finished")
