from __future__ import annotations
import logging
from mudata import read_h5mu, MuData
from sys import stdout
from collections.abc import Iterable

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

### VIASH START
par = {
    "input": "resources_test/concat_test_data/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix_subset.h5mu",
    "output": "foo.h5mu",
    "input_id": "mouse",
    "make_observation_keys_unique": True
}
### VIASH END


def make_observation_keys_unique(sample_id: str, sample: MuData) -> None:
    """
    Make the observation keys unique across all samples. At input,
    the observation keys are unique within a sample. By adding the sample name
    (unique for a sample) to each observation key, the observation key is made
    unique across all samples as well.
    """
    logger.info('Making observation keys unique across all samples.')
    sample.obs.index = f"{sample_id}_" + sample.obs.index
    make_observation_keys_unique_per_mod(sample_id, sample)


def make_observation_keys_unique_per_mod(sample_id: str, sample: MuData) -> None:
    """
    Updating MuData.obs_names is not allowed (it is read-only).
    So the observation keys for each modality has to be updated manually.
    """
    for mod in sample.mod.values():
        mod.obs_names = f"{sample_id}_" + mod.obs_names

def main():
    input_data = read_h5mu(par["input"])
    input_data.obs[par["obs_output"]] = par["input_id"]
    for mod_data in input_data.mod.values():
        mod_data.obs[par["obs_output"]] = par["input_id"]
    if par["make_observation_keys_unique"]:
        make_observation_keys_unique(par["input_id"], input_data)
    logger.info("Writing out data to '%s'.", par["output"])
    input_data.write_h5mu(par["output"], compression="gzip")

if __name__ == '__main__':
    main()