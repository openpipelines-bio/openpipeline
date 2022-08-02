import scanpy as sc
import muon as mu
import multiprocessing
import logging
from sys import stdout

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "output": "output.h5mu",
    "modality": ["rna"],
    "obs_keys": [],
}
meta = {"functionality_name": "lognorm"}
## VIASH END

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

logger.info("Reading input mudata")
mdata = mu.read_h5mu(par["input"])
mdata.var_names_make_unique()

if (
    par["obs_keys"] is not None
    and len(par["obs_keys"]) > 0
):

    for mod in par["modality"]:
        logger.info("Regress out variables on modality %s", mod)
        data = mdata.mod[mod]
        
        sc.pp.regress_out(
            data, 
            keys=par["obs_keys"], 
            n_jobs=multiprocessing.cpu_count() - 1
        )

# # can we assume execution_log exists?
# if mdata.uns is None or "execution_log" not in mdata.uns:
#     mdata.uns["execution_log"] = []
# # store new entry
# new_entry = {"component": meta["functionality_name"], "params": par}
# mdata.uns["execution_log"].append(new_entry)

logger.info("Writing to file")
mdata.write(filename=par["output"])
