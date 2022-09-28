import mudata
import scanpy as sc
import logging
from sys import stdout
import re
import pandas as pd
import numpy as np

# set logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

## VIASH START
par = {
  "sample_id": "foo", 
  "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5",
  "input_metrics_summary": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_metrics_summary.csv",
  "obs_sample_id": "sample_id",
  "obsm_metrics": "metrics_summary",
  "output": "foo.h5mu",
  "id_to_obs_names": True,
  "min_genes": 100,
  "min_counts": 1000
}
## VIASH END

logger.info("Reading %s.", par["input"])
adata = sc.read_10x_h5(par["input"], gex_only=False)

# store sample_id in .obs
if par["sample_id"] and par["obs_sample_id"]:
  logger.info("Storing sample_id '%s' in .obs['%s]'.", par['sample_id'], par['obs_sample_id'])
  adata.obs[par["obs_sample_id"]] = par["sample_id"]

# combine sample_id and barcode in obs_names
if par["sample_id"] and par["id_to_obs_names"]:
  logger.info("Combining obs_names and sample_id")
  # strip the number from '<10x_barcode>-<number>'
  replace = re.compile('-\\d+$')
  adata.obs_names = [ replace.sub('', obs_name) + "_" + par["sample_id"] for obs_name in adata.obs_names ]

# set the gene ids as var_names
logger.info("Renaming var columns")
adata.var = adata.var\
  .rename_axis("gene_symbol")\
  .reset_index()\
  .set_index("gene_ids")

# parse metrics summary file and store in .obsm or .obs
if par["input_metrics_summary"]:
  logger.info("Reading metrics summary file '%s'", par['input_metrics_summary'])

  def read_percentage(val):
      try:
          return float(val.strip('%')) / 100
      except AttributeError:
          return val

  df = pd.read_csv(par["input_metrics_summary"], decimal=".", quotechar='"', thousands=",").applymap(read_percentage)
  metrics_summary = df.iloc[np.repeat(0, adata.n_obs)]
  metrics_summary.index = adata.obs_names
    
  if par["obsm_metrics"]:
    logger.info("Storing metrics summary in .obs['%s']", par['obsm_metrics'])
    adata.obsm[par["obsm_metrics"]] = metrics_summary
  else:
    logger.info("Storing metrics summary in .obs")
    adata.obs = adata.obs.join(metrics_summary)

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

logger.info("Writing %s", par["output"])
mdata.write_h5mu(par["output"])
