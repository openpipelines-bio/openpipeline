import sys
from functools import partial
from operator import setitem

import numpy as np
from anndata import AnnData
from mudata import read_h5ad
from muon import prot as pt

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "modality": "prot",
    "output": "foo.h5mu",
    "layer": None,
}
## VIASH END


sys.path.append(meta["resources_dir"])
from compress_h5mu import write_h5ad_to_h5mu_with_compression


def main():
    input_data = read_h5ad(par["input"], mod=par["modality"])
    counts = (
        input_data.layers[par["input_layer"]] if par["input_layer"] else input_data.X
    )
    number_of_observations = input_data.n_obs
    assert number_of_observations > 1, (
        f"Need at least two observations to perform CLR normalization. Found: {number_of_observations}"
    )
    # muon's clr fails for integer counts, so cast to float first.
    # See https://github.com/scverse/muon/issues/257
    if not np.issubdtype(counts.dtype, np.floating):
        counts = counts.astype(np.float32)
    # CLR always normalizes the .X layer, so we have to create an AnnData object with
    # the input layer at .X. AnnData(X=counts) shares memory with input_data, so
    # inplace=False is required to leave the input layer untouched for the output.
    # The flavor is set explicitly so that a change to the muon default does not
    # silently alter the normalized counts.
    normalized_counts = pt.pp.clr(
        AnnData(X=counts), axis=par["axis"], inplace=False, flavor="seurat"
    )
    if not normalized_counts:
        raise RuntimeError("CLR failed to return the requested output layer")

    output_layer_setter = (
        partial(setattr, input_data, "X")
        if not par["output_layer"]
        else partial(setitem, input_data.layers, par["output_layer"])
    )
    output_layer_setter(normalized_counts.X)
    write_h5ad_to_h5mu_with_compression(
        par["output"],
        par["input"],
        par["modality"],
        input_data,
        par["output_compression"],
    )


if __name__ == "__main__":
    main()
