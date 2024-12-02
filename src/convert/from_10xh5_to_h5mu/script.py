import mudata
import scanpy as sc
import sys
import pandas as pd

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_raw_feature_bc_matrix.h5",
    "input_metrics_summary": None,
    "uns_metrics": "metrics_cellranger",
    "output": "foo.h5mu",
    "min_genes": None,
    "min_counts": None,
    "output_compression": "gzip",
}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
logger = setup_logger()

logger.info("Reading %s.", par["input"])
adata = sc.read_10x_h5(par["input"], gex_only=False)

# set the gene ids as var_names
logger.info("Renaming var columns")
adata.var = adata.var\
    .rename_axis("gene_symbol")\
    .reset_index()\
    .set_index("gene_ids")

# parse metrics summary file and store in .uns
if par["input_metrics_summary"] and par["uns_metrics"]:
    logger.info("Reading metrics summary file '%s'", par['input_metrics_summary'])

    def read_percentage(val):
        try:
            return float(val.strip('%')) / 100
        except AttributeError:
            return val

    metrics_summary = pd.read_csv(par["input_metrics_summary"], decimal=".", quotechar='"', thousands=",").applymap(read_percentage)

    logger.info("Storing metrics summary in .uns['%s']", par['uns_metrics'])
    adata.uns[par["uns_metrics"]] = metrics_summary
else:
    is_none = "input_metrics_summary" if not par["input_metrics_summary"] else "uns_metrics"
    logger.info("Not storing metrics summary because par['%s'] is None", is_none)

# might perform basic filtering to get rid of some data
# applicable when starting from the raw counts
if par["min_genes"]:
    logger.info("Filtering with min_genes=%d", par['min_genes'])
    sc.pp.filter_cells(adata, min_genes=par["min_genes"])

if par["min_counts"]:
    logger.info("Filtering with min_counts=%d", par['min_counts'])
    sc.pp.filter_cells(adata, min_counts=par["min_counts"])

# generate output
logger.info("Convert to mudata")
mdata = mudata.MuData(adata)

# override root .obs and .uns
mdata.obs = adata.obs
mdata.uns = adata.uns

# write output
logger.info("Writing %s", par["output"])
mdata.write_h5mu(par["output"], compression=par["output_compression"])
