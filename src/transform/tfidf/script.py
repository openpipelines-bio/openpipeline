import sys
import mudata
import muon

## VIASH START
par = {
    "input": "work/d9/3adbd080e0de618d44b59b1ec81685/run.output.h5mu",
    "output": "output.h5mu",
    "scale_factor": 10000,
    "modality": "atac",
    "input_layer": None,
    "output_layer": None,
    "output_compression": "gzip",
    "log_idf": True,
    "log_tf": True,
    "log_tfidf": False
}
meta = {"functionality_name": "lognorm"}
## VIASH END

sys.path.append(meta["resources_dir"])
# START TEMPORARY WORKAROUND setup_logger
# reason: resources aren't available when using Nextflow fusion
# from setup_logger import setup_logger
def setup_logger():
    import logging
    from sys import stdout

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    console_handler = logging.StreamHandler(stdout)
    logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
    console_handler.setFormatter(logFormatter)
    logger.addHandler(console_handler)

    return logger
# END TEMPORARY WORKAROUND setup_logger
logger = setup_logger()

logger.info("Reading input mudata")
mdata = mudata.read_h5mu(par["input"])
mdata.var_names_make_unique()

logger.info(par)

mod = par["modality"]
logger.info("Performing TF-IDF normalization on modality %s", mod)
adata = mdata.mod[mod]

muon.atac.pp.tfidf(
    adata,
    log_tf=par["log_tf"],
    log_idf=par["log_idf"],
    log_tfidf=par["log_tfidf"],
    scale_factor=par["scale_factor"],
    inplace=True,
    copy=False,
    from_layer=par["input_layer"],
    to_layer=par["output_layer"],
)

mdata.mod[mod] = adata

logger.info("Writing to file")
mdata.write_h5mu(filename=par["output"], compression=par["output_compression"])
