import sys
import os
import logging
import cellxgene_census
import mudata as mu
import anndata as ad
from scipy.sparse import csr_matrix
import obonet
import networkx

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
    "cell_query": "is_primary_data == True and cell_type_ontology_term_id in ['CL:0000136', 'CL:1000311', 'CL:0002616', 'CL:0002617', 'CL:0002615', 'CL:0001070', 'CL:1000310', 'CL:0000449', 'CL:1000309', 'CL:0002521', 'CL:0000448'] and assay_ontology_term_id in ['EFO:0009922', 'EFO:0009899', 'EFO:0011025', 'EFO:0030004', 'EFO:0030003'] and development_stage_ontology_term_id in ['HsapDv:0000084', 'HsapDv:0000098', 'HsapDv:0000097', 'HsapDv:0000099', 'HsapDv:0000096'] and suspension_type == 'cell'",
    "output": "output.h5mu",
    "output_compression": "gzip",
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


def write_mudata(mdata, output_location, compression):
    logger.info("Writing %s", output_location)
    


def main():

    census_connection = connect_census(
        par["input_database"],
        par["cellxgene_release"]
        ) 

    query_data = cellxgene_census.get_anndata(
        census = census_connection,
        obs_value_filter = par["cell_query"],
        organism = par["species"]
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
