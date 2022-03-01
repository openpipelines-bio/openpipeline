import muon
import scanpy as sc
import numpy as np

### VIASH START
par = {
    "input": "resources/test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "output": "/dev/null",
    "modality": ["rna"],
    "min_counts": int("200"),
    "max_counts": int("5000000"),
    "min_genes_per_cell": int("200"),
    "max_genes_per_cell": int("1500000"),
    "min_cells_per_gene": int("3"),
    "min_fraction_mito": float("0.0"),
    "max_fraction_mito": float("0.2"),
}
meta = {"functionality_name": "filter_on_umi_genes_and_mito"}
### VIASH END

print("Reading input data")
mudata = muon.read_h5mu(par["input"])

mudata.var_names_make_unique()

# sc.pp.filter_cells(data, min_genes_per_cell=0) # Hack to generate the n_genes column
# sc.pp.filter_genes(data, min_cells=0)

for mod in par["modality"]:
    print(f"Processing modality '{mod}'")
    data = mudata.mod[mod]

    print(f"  Unfiltered data: {data}")

    print("  Computing .obs['n_counts']")
    data.obs["n_counts"] = np.ravel(np.sum(data.X, axis=1))

    if par["min_cells_per_gene"] is not None:
        print(
            f"  Removing genes with non-zero values in <{par['min_cells_per_gene']} cells."
        )
        n_vars_before = data.n_vars
        sc.pp.filter_genes(data, min_cells=par["min_cells_per_gene"])
        if n_vars_before != data.n_vars:
            print(f"    Removed {n_vars_before - data.n_vars} genes.")

    print("  Removing cells with zero counts.")
    n_obs_before = data.n_obs
    sc.pp.filter_cells(data, min_genes=1)  # Hack to generate the n_genes column
    if n_obs_before != data.n_obs:
        print(f"    Removed {n_obs_before - data.n_obs} cells.")

    print(f"  Computing .obs['pct_mito']")
    mito_genes = data.var_names.str.contains("^[mM][tT]-")
    data.obs["pct_mito"] = np.ravel(
        np.sum(data[:, mito_genes].X, axis=1) / np.sum(data.X, axis=1)
    )

    if par["min_genes_per_cell"] is not None:
        print(f"  Removing cells with <{par['min_genes_per_cell']} expressed genes.")
        n_obs_before = data.n_obs
        sc.pp.filter_cells(data, min_genes=par["min_genes_per_cell"])
        if n_obs_before != data.n_obs:
            print(f"    Removed {n_obs_before - data.n_obs} cells.")

    if par["max_genes_per_cell"] is not None:
        print(f"  Removing cells with >{par['max_genes_per_cell']} expressed genes.")
        n_obs_before = data.n_obs
        sc.pp.filter_cells(data, max_genes=par["max_genes_per_cell"])
        if n_obs_before != data.n_obs:
            print(f"    Removed {n_obs_before - data.n_obs} cells.")

    if par["min_counts"] is not None:
        print(f"  Removing cells with <{par['min_counts']} total counts.")
        n_obs_before = data.n_obs
        data = data[data.obs["n_counts"] >= par["min_counts"], :]
        if n_obs_before != data.n_obs:
            print(f"    Removed {n_obs_before - data.n_obs} cells.")

    if par["max_counts"] is not None:
        print(f"  Removing cells with >{par['max_counts']} total counts.")
        n_obs_before = data.n_obs
        data = data[data.obs["n_counts"] <= par["max_counts"], :]
        if n_obs_before != data.n_obs:
            print(f"    Removed {n_obs_before - data.n_obs} cells.")

    if par["min_fraction_mito"] is not None:
        print(
            f"  Removing cells with <{par['min_fraction_mito']} percentage mitochondrial reads."
        )
        n_obs_before = data.n_obs
        data = data[data.obs["pct_mito"] <= par["max_fraction_mito"], :]
        if n_obs_before != data.n_obs:
            print(f"    Removed {n_obs_before - data.n_obs} cells.")

    if par["max_fraction_mito"] is not None:
        print(
            f"  Removing cells with >{par['max_fraction_mito']} percentage mitochondrial reads."
        )
        n_obs_before = data.n_obs
        data = data[data.obs["pct_mito"] >= par["min_fraction_mito"], :]
        if n_obs_before != data.n_obs:
            print(f"    Removed {n_obs_before - data.n_obs} cells.")

    data = data.copy()

    print("  Removing cells with zero counts (again).")
    n_obs_before = data.n_obs
    sc.pp.filter_genes(data, min_counts=1)
    if n_obs_before != data.n_obs:
        print(f"    Removed {n_obs_before - data.n_obs} cells.")

    print(f"  Filtered data: {data}")

# can we assume execution_log exists?
if mudata.uns is None or "execution_log" not in mudata.uns:
    mudata.uns["execution_log"] = []
# store new entry
new_entry = {"component": meta["functionality_name"], "params": par}
mudata.uns["execution_log"] = str(mudata.uns["execution_log"] + [new_entry])


print("Writing output data")
mudata.write(par["output"])
