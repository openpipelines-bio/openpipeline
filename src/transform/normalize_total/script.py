import scanpy as sc
import muon as mu
import logging 
from sys import stdout
## VIASH START
par = {
    "input": "work/d9/3adbd080e0de618d44b59b1ec81685/run.output.h5mu",
    "output": "output.h5mu",
    "target_sum": 10000,
    "modality": ["rna"],
    "exclude_highly_expressed": False
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

logger.info(par)

for mod in par["modality"]:
    logger.info("Performing total normalization on modality %s", mod)
    dat = mdata.mod[mod]
    logger.info(dat)
    sc.pp.normalize_total(
        dat,
        # target_sum=par["target_sum"],
        # exclude_highly_expressed=par["exclude_highly_expressed"]
    )

logger.info("Writing to file")
mdata.write(filename=par["output"])
