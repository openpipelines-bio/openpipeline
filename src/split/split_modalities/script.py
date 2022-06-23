from __future__ import annotations
import logging
import mudata as md
from sys import stdout
from pathlib import Path

### VIASH START
par = {
    "input": "./resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "output": "foo/",
    "compression": "gzip",
}
### VIASH END


logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)


def main() -> None:
    output_dir = Path(par["output"])
    if not output_dir.is_dir():
        output_dir.mkdir(parents=True)

    logger.info('Reading input file %s', par['input'])
    sample = md.read_h5mu(par["input"].strip())
    input_file = Path(par["input"])

    logger.info('Splitting up modalities %s', ", ".join(sample.mod.keys()))
    for mod_name, mod in sample.mod.items():
        new_sample = md.MuData({mod_name: mod})
        logger.info('Writing to %s_%s.h5mu', input_file.stem, mod_name)
        new_sample.write(output_dir / f"{input_file.stem}_{mod_name}.h5mu")

    logger.info("Finished")


if __name__ == "__main__":
    main()
