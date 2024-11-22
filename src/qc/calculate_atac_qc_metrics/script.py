import sys

import scanpy as sc
import muon as mu
from muon import atac as ac  # the module containing function for scATAC data processing
import numpy as np

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "fragments_path": None,
    "output": "foo.h5mu",
    "modality": "atac",
    "layer": None,
    "n_fragments_for_nucleosome_signal": 10e4,
    "n_tss": 3000,
    "nuc_signal_threshold": 2,
    "tss_threshold": 1.5,
}
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

def main():
    logger.info("Reading input data")
    mdata = mu.read(par["input"])

    atac = mdata.mod[par["modality"]]

    logger.info("Checking if QC columns are already calculated")
    for col in ("n_features_per_cell", "total_fragment_counts"):
        if col in atac.obs:
            logger.warning(f"{col} is already in atac.obs, dropping")
            atac.obs.drop(col, axis=1, inplace=True)
    
    logger.info("Calculating QC metrics")
    qc_metrics = sc.pp.calculate_qc_metrics(atac, percent_top=None, log1p=False, inplace=False, layer=par["layer"])

    logger.debug("Putting QC columns to ATAC adata")
    atac.obs["n_genes_by_counts"] = qc_metrics["n_features_per_cell"]
    atac.obs["total_counts"] = qc_metrics["total_fragment_counts"]

    logger.debug("Adding log-transformed total fragment counts")
    # log-transform total counts and add as column
    atac.obs["log_total_fragment_counts"] = np.log10(atac.obs["total_fragment_counts"])

    if par["fragments_path"] is not None:
        logger.info("Trying to locate frafments")
        ac.tl.locate_fragments(atac, fragments=par["fragments_path"])
    else:
        logger.info("Skipping fragment location: `fragments_path` is not set")
    
    # Calculate the nucleosome signal across cells
    if "files" in atac.uns and "fragments" in atac.uns["files"]:
        logger.info("Trying to calculate nucleosome signal")
        ac.tl.nucleosome_signal(atac, n=par["n_fragments_for_nucleosome_signal"] * atac.n_obs)
        atac.obs["nuc_signal_filter"] = [
            "NS_FAIL" if ns > par["nuc_signal_threshold"] else "NS_PASS"
            for ns in atac.obs["nucleosome_signal"]
        ]
    else:
        logger.info("Skipping nucleosome signal calculation: fragments information is not found")

    # If interval information is available, calculate TSS enrichment
    if "peak_annotation" in mdata.mod["atac"].uns.keys():
        tss = ac.tl.tss_enrichment(mdata, features=mdata.mod["atac"].uns["peak_annotation"],n_tss=par["n_tss"], random_state=666)
    
        tss.obs["tss_filter"] = [
            "TSS_FAIL" if score < par["tss_threshold"] else "TSS_PASS"
            for score in atac.obs["tss_score"]
        ]
    else:
        logger.info("Skipping TSS enrichment calculation: genes intervals are not found")

    logger.info("Writing output")
    mdata.write(par["output"], compression=par["output_compression"])
            
if __name__ == "__main__":
    main()