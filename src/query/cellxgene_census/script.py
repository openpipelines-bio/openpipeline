import sys
import os
import cellxgene_census
import mudata as mu
import anndata as ad

## VIASH START
par = {
    "input_database": "CellxGene",
    "modality": "rna",
    "cellxgene_release": "2023-05-15",
    "species": "homo_sapiens",
    "cell_query": "is_primary_data == True and cell_type_ontology_term_id in ['CL:0000136', 'CL:1000311', 'CL:0002616'] and suspension_type == 'cell'",
    "cells_filter_columns": ["dataset_id", "tissue", "assay", "disease", "cell_type"],
    "min_cells_filter_columns": 100,
    "output": "output.h5mu",
    "output_compression": "gzip",
}

meta = {'resources_dir': os.path.abspath('./src/query/cellxgene_census/')}
### VIASH END

sys.path.append(meta["resources_dir"])
# START TEMPORARY WORKAROUND setup_logger
# reason: resources aren't available when using Nextflow fusion
# from setup_logger import setup_logger
def setup_logger():
    import logging
    from sys import stdout

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    console_handler = logging.StreamHandler(stdout)
    logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
    console_handler.setFormatter(logFormatter)
    logger.addHandler(console_handler)

    return logger
# END TEMPORARY WORKAROUND setup_logger
logger = setup_logger()

def connect_census(input_database, release):
    """
    Connect to CellxGene Census or user-provided TileDBSoma object
    """
    if input_database != "CellxGene":
        raise NotImplementedError(
            "Custom census database is not implemented yet!"
            )

    logger.info(
        "Initializing %s release %s",
        input_database, release
        )
    return cellxgene_census.open_soma(
        census_version = release
        )


def get_anndata(census_connection, cell_query, species):
    logger.info(
        "Getting gene expression data based on %s query.",
        cell_query
        )
    return cellxgene_census.get_anndata(
        census = census_connection,
        obs_value_filter = cell_query,
        organism = species
    )


def add_cellcensus_metadata_obs(census_connection, query_data):
    logger.info(
    "Adding extented metadata to gene expression data."
    )
    census_datasets = census_connection["census_info"]["datasets"].read().concat().to_pandas()
    
    query_data.obs.dataset_id = query_data.obs.dataset_id.astype("category")

    dataset_info = census_datasets[census_datasets.dataset_id.isin(query_data.obs.dataset_id.cat.categories)]\
    [['collection_id', 'collection_name', 'collection_doi', 'dataset_id', 'dataset_title']]\
    .reset_index(drop=True)\
    .apply(lambda x: x.astype('category'))

    return query_data.obs.merge(
        dataset_info, on='dataset_id', how = 'left'
        )


def cellcensus_cell_filter(query_data, cells_filter_columns, min_cells_filter_columns):
    t0 = query_data.shape
    query_data = query_data[
        query_data.obs.groupby(cells_filter_columns)["soma_joinid"].transform('count') >= min_cells_filter_columns
        ]
    t1 = query_data.shape
    logger.info(
        'Removed %s cells based on %s min_cells_filter_columns of %s cells_filter_columns.'
        % ((t0[0] - t1[0]), min_cells_filter_columns, cells_filter_columns)
        )
    return query_data


def write_mudata(mdata, output_location, compression):
    logger.info("Writing %s", output_location)
    mdata.write_h5mu(
        output_location,
        compression=compression
        )


def main():

    # start dev
    logger.info('cells_filter_columns: %s' % par["cells_filter_columns"])
    logger.info('min_cells_filter_columns: %s' % par["min_cells_filter_columns"])
    # end dev
    
    census_connection = connect_census(
        par["input_database"],
        par["cellxgene_release"]
        ) 

    query_data = get_anndata(
        census_connection,
        par["cell_query"],
        par["species"]
        )
    
    query_data.obs = add_cellcensus_metadata_obs(
        census_connection,
        query_data
        )

    census_connection.close()
    del census_connection

    if par["cells_filter_columns"]:
        if not par["min_cells_filter_columns"]:
            raise NotImplementedError(
            "You specified cells_filter_columns, thus add min_cells_filter_columns!"
            )
        query_data = cellcensus_cell_filter(
            query_data,
            par["cells_filter_columns"],
            par["min_cells_filter_columns"]
            )
        
    query_data.var_names = query_data.var["feature_id"]
    query_data.var["gene_symbol"] = query_data.var["feature_name"]

    # Create empty mudata file
    mdata = mu.MuData({par["modality"]: ad.AnnData()})
    
    write_mudata(
        mdata,
        par["output"],
        par["output_compression"]
        )

    mu.write_h5ad(par["output"], data=query_data, mod=par["modality"])

if __name__ == "__main__":
    main()
