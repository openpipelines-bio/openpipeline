import liana
import mudata
# TODO: Remove when grouping labels exist
# For sign/PCA/
import numpy as np

### VIASH START
par = {
    "input": "integration.pca.output.h5mu",
    "output": "foo.h5mu",
    "compression": "gzip",
    "modality": "rna",
    "groupby": "bulk_labels",
    "resource_name": "consensus",
    "expr_prop": 0.1,
    "min_cells": 5,
    "aggregate_method": "rra",
    "return_all_lrs": False,
    "layer": "log_normalized",
    "n_perms": 100,
}
### VIASH END


def main():

    # Get input data
    mdata = mudata.read(par['input'].strip())
    mod = mdata.mod[par['modality']]

    # TODO: Remove when grouping labels exist
    # Add grouping labels
    foo = mod.obsm.to_df().iloc[:, 0]
    mod.obs[par['groupby']] = np.sign(foo).astype('category')

    # Solve gene labels
    orig_gene_label = mod.var.index
    mod.var_names = mod.var['gene_symbol'].astype(str)
    mod.var_names_make_unique()

    liana.mt.rank_aggregate(
        adata = mod,
        groupby = par['groupby'],
        resource_name = par["resource_name"],
        expr_prop = par["expr_prop"],
        min_cells = par["min_cells"],
        aggregate_method = par["aggregate_method"],
        return_all_lrs = par["return_all_lrs"],
        layer = par["layer"],
        n_perms = par["n_perms"],
        verbose = True,
        inplace = True,
        use_raw = False
    )

    # Return original gene labels
    mod.var_names = orig_gene_label

    # TODO: make sure compression is needed
    mdata.write_h5mu(par['output'].strip(), compression=par['compression'])

if __name__ == "__main__":
    main()