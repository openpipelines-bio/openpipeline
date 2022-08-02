import json
import scanpy as sc
import scipy
import muon as mu
import anndata
import logging
from sys import stdout

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.csv",
    "output": "output.h5ad",
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

logger.info("Converting.")
data.var_names_make_unique()
data.X = scipy.sparse.csr_matrix(data.X)
data.raw = data

muon = mu.MuData({"rna": data})

for key, value in json.loads(par["conversions_obsm"]).items():
    if key in data.obsm:
        muon.mod[value] = anndata.AnnData(data.obsm[key])
        del muon["rna"].obsm[key]

logger.info("Writing %s.", par["output"])
muon.write_h5mu(filename=par["output"])
logger.info("Finished")
