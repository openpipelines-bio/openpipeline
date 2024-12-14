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
    "log_tfidf": False,
}
meta = {"name": "tfidf"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

logger.info("Reading input mudata")
mdata = mudata.read_h5mu(par["input"])

logger.info(par)

mod = par["modality"]
logger.info("Performing TF-IDF normalization on modality %s", mod)
adata = mdata.mod[mod].copy()

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

mdata.mod[mod].layers[par["output_layer"]] = adata.layers[par["output_layer"]]

logger.info("Writing to file")
mdata.write_h5mu(filename=par["output"], compression=par["output_compression"])
