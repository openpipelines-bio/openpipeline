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

# where to find the obo files
cl_obo_folder = "/opt/cellxgene_census_cl_ontology/"

## VIASH START
par = {
    "input_database": "CellxGene",
    "modality": "rna",
    "cell_ontology_release": "2023-05-22",
    "cellxgene_release": "2023-05-15",
    "species": "homo_sapiens",
    "cell_type": ["mesothelial fibroblast"],
    "tissue": None,
    "technology": None,
    "suspension": None,
    "is_primary_data": True,
    "output": "output.h5mu",
    "output_compression": "gzip",
    "metadata_only": True
}
### VIASH END


def connect_census():
    """Connect to CellxGene Census or user-provided TileDBSoma object

    Returns:
        _type_: _description_
    """
    if par["input_database"] != "CellxGene":
        logger.info(
            "Custom census database is not implemented yet!"
            )
        sys.exit(1)

    logger.info(
        "Initializing %s release %s",
        par["input_database"],par["cellxgene_release"]
        )
    return cellxgene_census.open_soma(
        census_version = par["cellxgene_release"]
        )


def read_cell_ontology():
    """Reads Cell Type OBO Foundry file

    Returns:
        graph: cell type obo ontology
    """
    return obonet.read_obo(
        f"{cl_obo_folder}/cl.obo",
        encoding="utf-8"
    )


def get_cell_type_terms():
    """_summary_

    Returns:
        dict: _description_
        list: _description_
    """
    logging.info("Reading Cell Ontology OBO")
    co = read_cell_ontology()
    id_to_name = {id_: data.get("name") for id_, data in co.nodes(data=True)}
    id_to_name = {k: v for k, v in id_to_name.items() if v is not None}
    name_to_id = {v: k for k, v in id_to_name.items()}

    # TODO: can be more pythonic
    lower_hierarchy_cell_of_interest_map = {}
    cell_of_interest_terms = []
    for cell_type in par["cell_type"]:
        logging.info("Locating all child cell types of %s", cell_type)
        node = name_to_id[cell_type]
        lower_hierarchy_cell_of_interest_map.update(
            {parent:id_to_name[parent] for parent, child, key in co.in_edges(node, keys=True) if key == 'is_a'}
            )
        lower_hierarchy_cell_of_interest_map.update({node: cell_type})
        cell_of_interest_terms.extend(
            list(lower_hierarchy_cell_of_interest_map.keys())
            )

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


def build_census_query():
    _query = f'is_primary_data == {par["is_primary_data"]}'
    if par["cell_type"]:
        _, cell_of_interest_terms = get_cell_type_terms()
        _query = _query + f' and cell_type_ontology_term_id in {cell_of_interest_terms}'
    if par["tissue"]:
        _query = _query + f' and tissue in {par["tissue"]}'
    if par["technology"]:
        _query = _query + f' and assay in {par["technology"]}'
    if par["suspension"]:
        _query = _query + f' and suspension_type in {par["suspension"]}'
    logging.info("Census query: %s", _query)

    return _query

def extract_metadata(census_connection, query):
    logging.info("Extracting only metadata")
    return (
        census_connection["census_data"][par["species"]]
        .obs.read(
            value_filter=query
        )
        .concat()
        .to_pandas()
        )


def extract_metadata_expression(census_connection, query):
    logging.info("Extracting metadata and gene expression matrix")
    return (
        cellxgene_census.get_anndata(
            census_connection,
            organism = par["species"],
            obs_value_filter = query
            )
        )


def write_mudata(mdata):
    logger.info("Writing %s", par["output"])
    mdata.write_h5mu(
        par["output"],
        compression = par["output_compression"]
        )


def main():

    census_connection = connect_census()
    query = build_census_query()

    if par["metadata_only"]:
        query_data = extract_metadata(census_connection, query)

    else:
        query_data = extract_metadata_expression(census_connection, query)

    # clean-up mem
    census_connection.close()
    del census_connection

    if isinstance(query_data, pd.core.frame.DataFrame):
        mdata = mu.MuData(
            {par["modality"]: ad.AnnData(obs=query_data)}
            )
    else:
        query_data.X = csr_matrix(query_data.X)
        mdata = mu.MuData(
            {par["modality"]: query_data}
            )
        mdata.var_names_make_unique()

    write_mudata(mdata)

if __name__ == "__main__":
    main()
