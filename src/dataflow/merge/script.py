from __future__ import annotations
from pathlib import Path
from sys import stdout
import logging
import mudata as md



### VIASH START
par = {
    "input": ["./resources_test/merge/pbmc_1k_protein_v3_filtered_feature_bc_matrix_rna.h5mu",
              "./resources_test/merge/pbmc_1k_protein_v3_filtered_feature_bc_matrix_prot.h5mu"],
    "output": "foo.h5mu"

}
### VIASH END

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)


def main():
    logger.info('Reading input files %s', ",".join(par["input"]))
    input_samples = [md.read_h5mu(path) for path in par["input"]]

    logger.info('Merging into single object.')
    sample_modalities = {}
    for input_sample in input_samples:
        for mod_name, mod_data in input_sample.mod.items():
            if mod_name in sample_modalities:
                raise ValueError(f"Modality '{mod_name}' was found in more than 1 sample.")
            sample_modalities[mod_name] = mod_data

    merged = md.MuData(sample_modalities)
    merged.write_h5mu(par["output"], compression="gzip")
    logger.info('Finished')


if __name__ == '__main__':
    main()