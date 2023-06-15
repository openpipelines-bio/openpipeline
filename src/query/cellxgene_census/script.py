import sys
import logging
import cellxgene_census
import mudata as mu
import anndata as ad
import pandas as pd
from scipy.sparse import csr_matrix
import obonet

# set up logger
logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(sys.stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

## VIASH START
par = {
    "input_database": "CellxGene",
    "modality": "rna",
    "cellxgene_release": "2023-05-15",
    "species": "homo_sapiens",
    "cell_type": ["mesothelial fibroblast"],
    "tissue": None,
    "technology": None,
    "suspension": None,
    "is_primary_data": True,
    "obs_column_names": ["disease"],
    "output": "output.h5mu",
    "output_compression": "gzip",
    "metadata_only": True
}
### VIASH END


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


def read_cell_ontology():
    """Reads Cell Type OBO Foundry file

    Returns:
        graph: cell type obo ontology
    """
    return obonet.read_obo(
        "./cl-base.obo",
        encoding="utf-8"
    )


def get_cell_type_terms(cell_types):
    logger.info("Reading Cell Ontology OBO")

    co = read_cell_ontology()
    id_to_name = {
        id_: data.get("name")
        for id_, data in co.nodes(data=True)
        }
    id_to_name = {
        term_id: term_name
        for term_id, term_name in id_to_name.items()
        if term_name is not None
        }
    name_to_id = {
        term_name: term_id
        for term_id, term_name in id_to_name.items()
        }

    # TODO: can be more pythonic
    lower_hierarchy_cell_of_interest_map = {}
    cell_of_interest_terms = []
    for cell_type in cell_types:
        logger.info("Locating all child cell types of %s", cell_type)
        node = name_to_id[cell_type]
        lower_hierarchy_cell_of_interest_map.update(
            {
                parent:id_to_name[parent]
                for parent, child, key in co.in_edges(node, keys=True)
                if key == 'is_a'
                }
            )
        lower_hierarchy_cell_of_interest_map.update({node: cell_type})
        cell_of_interest_terms.extend(
            list(lower_hierarchy_cell_of_interest_map.keys())
            )

    logger.info(lower_hierarchy_cell_of_interest_map)

    return lower_hierarchy_cell_of_interest_map, cell_of_interest_terms


# TODO: function to explore cell types available in query data
# def view_available_cell_types(lower_hierarchy_cell_of_interest_map, cell_of_interest_terms):
#     cells_of_interest_query_terms = list(adata.obs.cell_type_ontology_term_id.unique())
#     avail_cells_overview = []
#     for term in cell_of_interest_terms:
#         if term in cells_of_interest_query_terms:
#             avail_cells_overview.append("Available: {}-{}".format(term, lower_hierarchy_cell_of_interest_map[term]))
#         else:
#             avail_cells_overview.append("Unavailable: {}-{}".format(term, lower_hierarchy_cell_of_interest_map[term]))
#     avail_cells_overview.sort()
#     logging.info(avail_cells_overview)


def build_census_query(par):
    _query = f'is_primary_data == {par["is_primary_data"]}'
    query_builder = {
        'cell_type': f' and cell_type_ontology_term_id in {get_cell_type_terms(par["cell_type"])}',
        'tissue': f' and tissue in {par["tissue"]}',
        'technology': f' and assay in {par["technology"]}',
        'suspension': f' and suspension_type in {par["suspension"]}',
    }
    for parameter_name, query_part in query_builder.items():
        if par[parameter_name]:
            _query += query_part

    return _query


def extract_metadata(census_connection, query, species, obs_column_names):
    logger.info("Extracting only metadata")

    query_data = census_connection["census_data"][species].obs.read(
        value_filter=query,
        column_names=obs_column_names).concat().to_pandas()

    return ad.AnnData(obs=query_data)


def extract_metadata_expression(
    census_connection, query,species, obs_column_names):
    logger.info("Extracting metadata and gene expression matrix")

    query_data = cellxgene_census.get_anndata(
        census_connection,
        organism = species,
        obs_value_filter = query,
        column_names = {"obs": obs_column_names}
        )

    query_data.X = csr_matrix(query_data.X)
    query_data.var_names_make_unique()

    return query_data


def write_mudata(mdata, output_location, compression):
    logger.info("Writing %s", output_location)
    mdata.write_h5mu(
        output_location,
        compression = compression
        )


def main():

    census_connection = connect_census(
        par["input_database"],
        par["cellxgene_release"]
        )
    query = build_census_query(par)

    if par["metadata_only"]:
        query_data = extract_metadata(
            census_connection, query, par["species"], par["obs_column_names"]
            )

    else:
        query_data = extract_metadata_expression(
            census_connection, query, par["species"], par["obs_column_names"]
            )

    census_connection.close()
    del census_connection

    mdata = mu.MuData(
        {par["modality"]: query_data}
        )

    write_mudata(mdata, par["output"], par["compression"])

if __name__ == "__main__":
    main()
