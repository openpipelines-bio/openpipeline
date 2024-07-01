from __future__ import annotations
import sys
import mudata as md
from sys import stdout
from pathlib import Path
import pandas as pd

### VIASH START
par = {
    "input": "./resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "output": "foo/",
    "output_types": "foo_types.csv",
    "output_compression": "gzip",
}
meta = {
    "resources_dir": "."
}
### VIASH END

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

def main() -> None:
    output_dir = Path(par["output"])
    logger.info("Creating output directory '%s' if it does not exist", output_dir)
    if not output_dir.is_dir():
        logger.info("Creating %s", output_dir)
        output_dir.mkdir(parents=True)

    logger.info("Reading input file '%s'", par['input'])
    input_file = Path(par["input"].strip())
    sample = md.read_h5mu(input_file)

    logger.info('Creating output types CSV.')
    modalities = list(sample.mod.keys())

    logger.info("Found the following modalities:\n%s", "\n".join(modalities))
    names = {mod_name: f"{input_file.stem}_{mod_name}.h5mu"
             for mod_name in modalities}
    output_files = list(names.values())
    logger.info("Will be creating the following output .h5mu files:\n%s", "\n".join(output_files))
    df = pd.DataFrame({"name": modalities, "filename": output_files})
    logger.info("Writing output_types CSV file to '%s'.", par["output_types"])
    df.to_csv(par["output_types"], index=False)

    logger.info('Splitting input file into unimodal output files.')
    for mod_name, mod in sample.mod.items():
        logger.info("Processing modality '%s'", mod_name)
        new_sample = md.MuData({mod_name: mod})
        logger.info("Writing to '%s', with compression '%s'", names[mod_name], par["output_compression"])
        new_sample.write_h5mu(output_dir / names[mod_name], compression=par["output_compression"])
        logger.info("Done writing output file.")
    logger.info("Finished")

if __name__ == "__main__":
    main()
