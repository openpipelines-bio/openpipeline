import sys
from muon import prot as pt
from mudata import read_h5ad
from anndata import AnnData
from functools import partial
from operator import setitem

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
    if par["input_layer"]:
        input_data = AnnData(X=input_data.layers[par["input_layer"]])
    # CLR always normalizes the .X layer, so we have to create an AnnData file with
    # the input layer at .X
    normalized_counts = pt.pp.clr(input_data, axis=par["axis"], inplace=False)
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
