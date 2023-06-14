from __future__ import annotations
import logging
import mudata as md
from sys import stdout
from pathlib import Path
import pandas as pd

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

### VIASH START
par = {
    "input": "./resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "output": "foo/",
    "output_types": "foo_types.csv",
    "compression": "gzip",
}
### VIASH END

def main() -> None:
    output_dir = Path(par["output"])
    if not output_dir.is_dir():
        output_dir.mkdir(parents=True)

    logger.info('Reading input file %s', par['input'])
    sample = md.read_h5mu(par["input"].strip())
    input_file = Path(par["input"])

    logger.info('Creating output types csv')

    names = {mod_name: f"{input_file.stem}_{mod_name}.h5mu"
        for mod_name in sample.mod.keys() }
    df = pd.DataFrame({"name": list(names.keys()), "filename": list(names.values())})
    df.to_csv(par["output_types"], index=False)

    logger.info('Splitting up modalities %s', ", ".join(sample.mod.keys()))
    for mod_name, mod in sample.mod.items():
        new_sample = md.MuData({mod_name: mod})
        logger.info('Writing to %s', names[mod_name])
        new_sample.write_h5mu(output_dir / names[mod_name], compression=par["output_compression"])

    logger.info("Finished")


if __name__ == "__main__":
    main()
