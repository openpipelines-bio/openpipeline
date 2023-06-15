from muon import prot as pt
from mudata import read_h5mu

## VIASH START
par = {
    'input': 'resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu',
    'modality': 'prot',
    'output': "foo.h5mu"
}
## VIASH END


def main():
    input_h5mu = read_h5mu(par['input'])
    modality = input_h5mu[par['modality']]
    normalized_counts = pt.pp.clr(modality, inplace=False if par['output_layer'] else True)
    if normalized_counts:
        input_h5mu[par["modality"]].layers[par['output_layer']] = normalized_counts.X
    input_h5mu.write_h5mu(par['output'], compression=par["output_compression"])

if __name__ == "__main__":
    main()