from __future__ import annotations
import sys
from mudata import read_h5mu, MuData

### VIASH START
par = {
    "input": "resources_test/concat_test_data/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix_subset.h5mu",
    "output": "foo.h5mu",
    "input_id": "mouse",
    "make_observation_keys_unique": True
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

def make_observation_keys_unique(sample_id: str, sample: MuData) -> None:
    """
    Make the observation keys unique across all samples. At input,
    the observation keys are unique within a sample. By adding the sample name
    (unique for a sample) to each observation key, the observation key is made
    unique across all samples as well.
    """
    logger.info("Making observation keys unique across all "
                "samples by appending prefix '%s' to the observation names.",
                sample_id)
    sample.obs.index = f"{sample_id}_" + sample.obs.index
    make_observation_keys_unique_per_mod(sample_id, sample)
    logger.info("Done making observation keys unique.")


def make_observation_keys_unique_per_mod(sample_id: str, sample: MuData) -> None:
    """
    Updating MuData.obs_names is not allowed (it is read-only).
    So the observation keys for each modality has to be updated manually.
    """
    for mod_name, mod in sample.mod.items():
        logger.info("Processing modality '%s'", mod_name)
        mod.obs_names = f"{sample_id}_" + mod.obs_names

def main():
    logger.info("Reading input file '%s'.", par["input"])
    input_data = read_h5mu(par["input"])
    logger.info("Adding column '%s' to global .obs dataframe, populated with ID '%s'", 
                par["obs_output"], par["input_id"])
    input_data.obs[par["obs_output"]] = par["input_id"]
    logger.info("Done adding column to global .obs")
    for mod_name, mod_data in input_data.mod.items():
        logger.info("Adding column '%s' to .obs dataframe for modality '%s', "
                    "populated with ID '%s'", par["obs_output"], mod_name, par["input_id"])
        mod_data.obs[par["obs_output"]] = par["input_id"]
    logger.info("Done adding per-modality columns.")
    if par["make_observation_keys_unique"]:
        make_observation_keys_unique(par["input_id"], input_data)
    logger.info("Writing out data to '%s' with compression '%s'.", 
                par["output"], par["output_compression"])
    input_data.write_h5mu(par["output"], compression=par["output_compression"])
    logger.info("Finished")

if __name__ == '__main__':
    main()