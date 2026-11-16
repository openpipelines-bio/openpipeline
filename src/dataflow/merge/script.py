from __future__ import annotations
import sys
import mudata as md


### VIASH START
par = {
    "input": [
        "./resources_test/merge/pbmc_1k_protein_v3_filtered_feature_bc_matrix_rna.h5mu",
        "./resources_test/merge/pbmc_1k_protein_v3_filtered_feature_bc_matrix_prot.h5mu",
    ],
    "output": "foo.h5mu",
    "output_compression": None,
}
meta = {"resources_dir": "src/utils"}
### VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()


def main():
    logger.info("Reading input files %s", ",".join(par["input"]))
    input_samples = [md.read_h5mu(path) for path in par["input"]]

    logger.info("Merging into single object.")
    sample_modalities = {}
    for input_sample in input_samples:
        for mod_name, mod_data in input_sample.mod.items():
            if mod_name in sample_modalities:
                raise ValueError(
                    f"Modality '{mod_name}' was found in more than 1 sample."
                )
            sample_modalities[mod_name] = mod_data

    merged = md.MuData(sample_modalities)
    merged.write_h5mu(par["output"], compression=par["output_compression"])
    logger.info("Finished")


if __name__ == "__main__":
    main()
