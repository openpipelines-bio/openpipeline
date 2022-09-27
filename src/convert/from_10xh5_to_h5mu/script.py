import mudata
import scanpy as sc
import logging
from sys import stdout
import re
import pandas as pd
import csv

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

if par["sample_id"] is not None and par["obs_sample_id"] is not None:
  logger.info(f"Storing sample_id '{par['sample_id']}' in .obs['{par['obs_sample_id']}]'.")
  adata.obs[par["obs_sample_id"]] = par["sample_id"]

if par["sample_id"] is not None and par["id_to_obs_names"] == True:
  logger.info("Combining obs_names and sample_id")
  replace = re.compile('-\\d+$')
  adata.obs_names = [ replace.sub('', obs_name) + "_" + par["sample_id"] for obs_name in adata.obs_names ]
  
logger.info("Renaming var columns")
adata.var = adata.var\
  .rename_axis("gene_symbol")\
  .reset_index()\
  .set_index("gene_ids")
  #.rename(columns={"gene_ids":"gene_id", "feature_types":"feature_type"})\

if par["input_metrics_summary"] is not None:
  logger.info(f"Reading metrics summary file '{par['input_metrics_summary']}'")
  # pd.read_csv(par["input_metrics_summary"], thousands=",")
  with open(par["input_metrics_summary"], newline='') as f:
    reader = csv.reader(f)
    data = list(reader)
    header = data[0]
    # header = [ s.lower().replace(' ', '_') for s in data[0] ]
    # ^ do we want to remove spaces in the names?
    
    def string2floatorint(x):
      y = re.sub("[^\\d\\.]", "", x)
      if '%' in x:
        return float(y) / 100
      elif '.' in x:
        return float(y)
      else:
        return int(y)
    
    values = [ string2floatorint(x) for x in data[1] ]
    
    metrics_summary = pd.DataFrame(dict(zip(header, values)), index=adata.obs_names)
    
  if par["obsm_metrics"] is not None:
    logger.info(f"Storing metrics summary in .obs['{par['obsm_metrics']}']")
    adata.obsm[par["obsm_metrics"]] = metrics_summary
  else:
    logger.info("Storing metrics summary in .obs")
    adata.obs = adata.obs.join(metrics_summary)

# might perform basic filtering to get rid of some data
# applicable when starting from the raw counts
if par["min_genes"] is not None:
  logger.info(f"Filtering with min_genes={par['min_genes']}")
  sc.pp.filter_cells(adata, min_genes=par["min_genes"])

if par["min_counts"] is not None:
  logger.info(f"Filtering with min_counts={par['min_counts']}")
  sc.pp.filter_cells(adata, min_counts=par["min_counts"])

logger.info("Convert to mudata")
mdata = mudata.MuData(adata)

logger.info("Writing %s.", par["output"])
mdata.write_h5mu(par["output"])
