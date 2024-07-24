import mudata as mu

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms_w_sample_id.h5mu",
    # "input": "resources_test/annotation_test_data/TS_Blood_filtered.h5mu",
    "modality": "rna",
    "var_key_input": "gene_symbol",
    # "var_key_input": None,
    "var_key_output": "genes",
    "output": "output.h5mu",
    "output_compression": "gzip"
}
## VIASH END

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

mdata = mu.read(par["input"])

if par["var_key_input"]:
    mdata.mod[par["modality"]].var.rename(columns={par["var_key_input"]: par["var_key_output"]}, inplace=True)
else:
    mdata.mod[par["modality"]].var[par["var_key_output"]] = mdata.mod[par["modality"]].var.index

# if par["var_key_input"]:
#     assert par["var_key_input"] not in mdata.mod[par["modality"]].var.columns
#     assert par["var_key_output"] in mdata.mod[par["modality"]].var.columns
# else:
#     assert par["var_key_output"] in mdata.mod[par["modality"]].var.columns

mdata.write_h5mu(par["output"], compression=par["output_compression"])
