import sys
import mudata as mu
import anndata as ad
import re

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "modality": "rna",
    "reference": "resources_test/annotation_test_data/TS_Blood_filtered.h5mu",
    "reference_obs_label": "gene_symbol",
    "output_query": "query_overlap.h5mu",
    "output_reference": "reference_overlap.h5mu",
}
meta = {
    "resources_dir": "src/labels_transfer/utils"
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

input_mudata = mu.read_h5mu(par["input"])
reference_mudata = mu.read_h5mu(par["reference"])

input_modality = input_mudata.mod[par["modality"]]
reference_modality = reference_mudata.mod[par["modality"]]

if par["var_query_gene_names"]:
    input_modality.var.index = [re.sub("\\.[0-9]+$", "", s) for s in input_modality.var[par["var_query_gene_names"]]]
    
if par["var_reference_gene_names"]:
    reference_modality.var.index = [re.sub("\\.[0-9]+$", "", s) for s in reference_modality.var[par["var_reference_gene_names"]]]

common_ens_ids = list(set(reference_modality.var.index).intersection(set(input_modality.var.index)))

reference_modality = reference_modality[:, common_ens_ids].copy()
input_modality = input_modality[:, common_ens_ids].copy()

input_mudata.mod[par["modality"]] = input_modality
reference_mudata.mod[par["modality"]] = reference_modality

input_mudata.write_h5mu(par["output_query"], compression=par["compression"])
reference_mudata.write_h5mu(par["output_reference"], compression=par["compression"])