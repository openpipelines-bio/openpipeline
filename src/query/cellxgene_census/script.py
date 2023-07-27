import sys
import os
import logging
import cellxgene_census
import mudata as mu
import anndata as ad
from scipy.sparse import csr_matrix
import obonet
import networx

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
    "metadata_only": False
}

meta = {'resources_dir': os.path.abspath('./src/query/cellxgene_census/')}

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


def read_cell_ontology(obo_file):
    """Reads Cell Type OBO Foundry file

    Returns:
        graph: cell type obo ontology
    """
    return obonet.read_obo(
        obo_file,
        encoding="utf-8"
    )


def get_child_terms_from_ontology(cell_types, ontology, relations_to_consider = None):
    
    if relations_to_consider:
        def filter_edge(n1, n2, e):
            return e in relations_to_consider
        
        ontology = networkx.subgraph_view(ontology, filter_edge = filter_edge)
    
    id_to_name = {id_: data.get('name') for id_, data in ontology.nodes(data=True)}
    name_to_id = {data['name']: id_ for id_, data in ontology.nodes(data=True) if 'name' in data}

    def flatten_extend(matrix):
        flat_list = []
        for row in matrix:
            flat_list.extend(row)
        return flat_list
    
    return set([cell_types[0]] + flatten_extend([sorted(subterm for subterm in networkx.ancestors(ontology, name_to_id[ct])) for ct in cell_types ]))

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


def build_census_query(par, obo_file):
    _query = f'is_primary_data == {par["is_primary_data"]}'
    query_builder = {
        'cell_type': f' and cell_type_ontology_term_id in {get_child_terms_from_ontology(cell_types = par["cell_type"], ontology = obo_file, relations_to_consider = ["is_a"])}',
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
    query = build_census_query(par, f"{meta['resources_dir']}/cl-base.obo")

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
    
    mdata["rna"].var_names = mdata["rna"].var["feature_id"]
    mdata["rna"].var["gene_symbol"] = mdata["rna"].var["feature_name"]

    write_mudata(mdata, par["output"], par["output_compression"])

if __name__ == "__main__":
    main()
