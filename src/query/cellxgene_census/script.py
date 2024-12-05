import sys
import cellxgene_census
import scanpy as sc
import mudata as mu

## VIASH START
par = {
    "input_uri": None,
    "census_version": "stable",
    "species": "mus_musculus",
    "obs_value_filter": "dataset_id == '49e4ffcc-5444-406d-bdee-577127404ba8'",
    "cell_filter_grouping": None,
    "cell_filter_minimum_count": None,
    "output": "output.h5ad",
    "output_modality": "rna",
    "output_compression": "gzip",
}
meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])


from setup_logger import setup_logger

logger = setup_logger()


def connect_census(uri, census_version):
    """
    Connect to CellxGene Census or user-provided TileDBSoma object
    """
    ver = census_version or "stable"
    logger.info(
        "Connecting to CellxGene Census at %s",
        f"'{uri}'" if uri else f"version '{ver}'",
    )
    return cellxgene_census.open_soma(uri=uri, census_version=ver)


def get_anndata(census_connection, par):
    logger.info(
        "Getting gene expression data based on `%s` query.", par["obs_value_filter"]
    )
    return cellxgene_census.get_anndata(
        census=census_connection,
        obs_value_filter=par["obs_value_filter"],
        organism=par["species"],
    )


def add_cellcensus_metadata_obs(census_connection, adata):
    logger.info("Adding additional metadata to gene expression data.")
    census_datasets = (
        census_connection["census_info"]["datasets"].read().concat().to_pandas()
    )

    adata.obs.dataset_id = adata.obs.dataset_id.astype("category")

    dataset_info = (
        census_datasets[
            census_datasets.dataset_id.isin(adata.obs.dataset_id.cat.categories)
        ][
            [
                "collection_id",
                "collection_name",
                "collection_doi",
                "dataset_id",
                "dataset_title",
            ]
        ]
        .reset_index(drop=True)
        .astype("category")
    )

    adata.obs = adata.obs.merge(dataset_info, on="dataset_id", how="left")


def filter_min_cells_per_group(adata, par):
    n_cells_before, _ = adata.shape
    cell_count = adata.obs.groupby(par["cell_filter_grouping"])[
        "soma_joinid"
    ].transform("count")
    adata = adata[cell_count >= par["cell_filter_minimum_count"]]
    n_cells_after, _ = adata.shape
    logger.info(
        "Removed %s cells based on %s cell_filter_minimum_count of %s cell_filter_grouping."
        % (
            (n_cells_before - n_cells_after),
            par["cell_filter_minimum_count"],
            par["cell_filter_grouping"],
        )
    )
    return adata


def filter_by_counts(adata, par):
    logger.info("Remove cells with few counts and genes with few counts.")
    n_cells_before, n_genes_before = adata.shape
    # remove cells with few counts and genes with few counts
    scanpy_proc = {
        par["cell_filter_min_counts"]: (sc.pp.filter_cells, "min_counts"),
        par["cell_filter_min_genes"]: (sc.pp.filter_cells, "min_genes"),
        par["gene_filter_min_counts"]: (sc.pp.filter_genes, "min_counts"),
        par["gene_filter_min_cells"]: (sc.pp.filter_genes, "min_cells"),
    }
    for threshold, (func, arg) in scanpy_proc.items():
        if threshold:
            func(adata, **{arg: threshold})
    n_cells_after, n_genes_after = adata.shape
    logger.info(
        "Removed %s cells and %s genes.",
        (n_cells_before - n_cells_after),
        (n_genes_before - n_genes_after),
    )


def move_x_to_layers(adata):
    logger.info("Move .X to .layers['counts']")
    adata.layers["counts"] = adata.X
    adata.X = None


def print_unique(adata, column):
    unique_values = adata.obs[column].unique().astype(str)
    formatted = "', '".join(unique_values[:50])
    if len(unique_values) > 50:
        formatted += ", ..."
    logger.info(f"Unique {column}: ['{formatted}']")


def print_summary(adata):
    logger.info(f"Resulting dataset: {adata}")

    logger.info("Summary of dataset:")
    for field in adata.obs.columns:
        print_unique(adata, field)


def write_anndata(adata, par):
    logger.info("Writing MuData object to '%s'", par["output"])

    mdata = mu.MuData({par["output_modality"]: adata})

    mdata.write_h5mu(par["output"], compression=par["output_compression"])


def main(par, meta):
    # check arguments
    if (par["cell_filter_grouping"] is None) != (
        par["cell_filter_minimum_count"] is None
    ):
        raise NotImplementedError(
            "You need to specify either both or none of the following parameters: cell_filter_grouping, cell_filter_minimum_count"
        )

    with connect_census(
        uri=par["input_uri"], census_version=par["census_version"]
    ) as conn:
        adata = get_anndata(conn, par)

        if par["add_dataset_metadata"]:
            add_cellcensus_metadata_obs(conn, adata)

    print(f"AnnData: {adata}", flush=True)

    if par["cell_filter_grouping"] is not None:
        adata = filter_min_cells_per_group(adata, par)

    # remove cells with few counts and genes with few counts
    filter_by_counts(adata, par)

    # logger.log(f"Filtered AnnData: {adata}")
    print(f"Filtered AnnData: {adata}", flush=True)

    # use feature_id as var_names
    adata.var_names = adata.var["feature_id"]

    # move .X to .layers["counts"]
    move_x_to_layers(adata)

    # print summary
    print_summary(adata)

    # write output to file
    write_anndata(adata, par)


if __name__ == "__main__":
    main(par, meta)
