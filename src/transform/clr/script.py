from muon import prot as pt
from mudata import read_h5mu
from anndata import AnnData
from functools import partial
from operator import setitem

## VIASH START
par = {
    'input': 'resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu',
    'modality': 'prot',
    'output': "foo.h5mu",
    'layer': None,
}
## VIASH END


def main():
    input_h5mu = read_h5mu(par['input'])
    modality = input_h5mu[par['modality']]
    input_data = modality
    if par["input_layer"]:
        input_data = AnnData(X=input_data.layers[par["input_layer"]])
    # CLR always normalizes the .X layer, so we have to create an AnnData file with
    # the input layer at .X
    normalized_counts = pt.pp.clr(input_data, axis=1, inplace=False)
    if not normalized_counts:
        raise RuntimeError("CLR failed to return the requested output layer")

    output_layer_setter = partial(setattr, modality, "X") \
                          if not par["output_layer"] \
                          else partial(setitem, modality.layers, par["output_layer"])
    output_layer_setter(normalized_counts.X)
    input_h5mu.write_h5mu(par['output'], compression=par["output_compression"])

if __name__ == "__main__":
    main()