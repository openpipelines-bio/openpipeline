from mudata import read_h5mu
from sys import stdout
import logging

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "output": "output.h5mu",
    "modality": ["rna"],
    "layer": ['test'],
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
    logger.info('Processing modalities: %s.', ','.join(par['modality']))
    for mod_name in par['modality']:
        mod = input_data.mod[mod_name]
        for layer in par['layer']:
            logger.info('Deleting layer %s from modality %s.', layer, mod_name)
            del mod.layers[layer]
    logger.info('Writing output to %s.', par['output'])
    input_data.write_h5mu(par['output'])
    logger.info('Finished.')

if __name__ == "__main__":
    main()