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
from compress_h5mu import write_h5ad_to_h5mu_with_compression

logger = setup_logger()

logger.info("Reading input modality %s from %s", par["modality"], par["input"])
adata = mudata.read_h5ad(par["input"], mod=par["modality"])

logger.info(par)

logger.info("Performing TF-IDF normalization")
input_data = adata.copy()

muon.atac.pp.tfidf(
    input_data,
    log_tf=par["log_tf"],
    log_idf=par["log_idf"],
    log_tfidf=par["log_tfidf"],
    scale_factor=par["scale_factor"],
    inplace=True,
    copy=False,
    from_layer=par["input_layer"],
    to_layer=par["output_layer"],
)

adata.layers[par["output_layer"]] = input_data.layers[par["output_layer"]]

logger.info("Writing to file")
write_h5ad_to_h5mu_with_compression(
    par["output"], par["input"], par["modality"], adata, par["output_compression"]
)
