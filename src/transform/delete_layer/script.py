from mudata import read_h5mu
from sys import stdout
import logging

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "layer": ['test'],
    "missing_ok": False
}
meta = {"functionality_name": "lognorm"}
## VIASH END

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

def main():
    logger.info('Reading input file %s.', par['input'])
    input_data = read_h5mu(par['input'])
    logger.info('Processing modality: %s.', par['modality'])
    mod_name = par['modality']
    mod = input_data.mod[mod_name]
    for layer in par['layer']:
        if layer not in mod.layers:
            if par['missing_ok']:
                continue
            raise ValueError(f"Layer '{layer}' is not present in modality {mod_name}.")
        logger.info('Deleting layer %s from modality %s.', layer, mod_name)
        del mod.layers[layer]
    logger.info('Writing output to %s.', par['output'])
    input_data.write_h5mu(par['output'], compression=par["output_compression"])
    logger.info('Finished.')

if __name__ == "__main__":
    main()