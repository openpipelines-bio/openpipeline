import sys
import warnings

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
    mdata = mu.read(par["input"])

    mdata.var_names_make_unique()
    atac = mdata.mod[par["modality"]]

    for col in ("n_features_per_cell", "total_fragment_counts"):
        if col in atac.obs:
            warnings.warn(f"{col} is already in atac.obs, dropping")
            atac.obs.drop(col, axis=1, inplace=True)
    
    # Calculate general qc metrics using scanpy
    sc.pp.calculate_qc_metrics(atac, percent_top=None, log1p=False, inplace=True, layer=par["layer"])

    # Rename columns
    atac.obs.rename(
        columns={
            "n_genes_by_counts": "n_features_per_cell",
            "total_counts": "total_fragment_counts",
        },
        inplace=True,
    )

    # log-transform total counts and add as column
    atac.obs["log_total_fragment_counts"] = np.log10(atac.obs["total_fragment_counts"])

    if par["fragments_path"] is not None:
        ac.tl.locate_fragments(atac, fragments=par["fragments_path"])
    
    # Calculate the nucleosome signal across cells
    if "files" in atac.uns and "fragments" in atac.uns["files"]:
        ac.tl.nucleosome_signal(atac, n=par["n_fragments_for_nucleosome_signal"] * atac.n_obs)
        atac.obs["nuc_signal_filter"] = [
            "NS_FAIL" if ns > par["nuc_signal_threshold"] else "NS_PASS"
            for ns in atac.obs["nucleosome_signal"]
        ]


    # If interval information is available, calculate TSS enrichment
    if "interval" in mdata.var or ("rna" in mdata.mod and "interval" in mdata.mod["rna"]):
        tss = ac.tl.tss_enrichment(mdata, n_tss=par["n_tss"], random_state=666)
    
        tss.obs["tss_filter"] = [
            "TSS_FAIL" if score < par["tss_threshold"] else "TSS_PASS"
            for score in atac.obs["tss_score"]
        ]

    mdata.write(par["output"], compression=par["output_compression"])
            
if __name__ == "__main__":
    main()