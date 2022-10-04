import mudata as mu
# import scanpy as sc
import logging
import tempfile
import subprocess
import os
import sys

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(sys.stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "total_droplets_included": 50000,
    "epochs": 150,
    "fpr": 0.01,
    "exclude_antibody_capture": False,
    "learning_rate": 0.001,
    "layer_corrected": "corrected",
    "cuda": False
}
meta = { 'temp_dir': 'foo', 'resources_dir': 'src/correction/cellbender_remove_background' }
## VIASH END

sys.path.append(meta['resources_dir'])
from helper import anndata_from_h5

logger.info("Reading input mudata")
mdata = mu.read_h5mu(par["input"])

mod = par["modality"]
logger.info("Performing log transformation on modality %s", mod)
data = mdata.mod[mod]

with tempfile.NamedTemporaryFile(suffix = ".h5ad") as input_file:
  with tempfile.TemporaryDirectory(prefix="cellbender-", dir=meta["temp_dir"]) as temp_dir:
    data.write_h5ad(input_file.name)
    
    output_file = os.path.join(temp_dir, "output.h5ad")
    output_report = os.path.join(temp_dir, "output.pdf")
    output_cell_barcodes = os.path.join(temp_dir, "output_cell_barcodes.csv")
    cmd_pars = [
      "cellbender",
      "remove-background",
      "--input", input_file.name,
      "--output", output_file
    ]

    if par["model"]:
      cmd_pars += [ "--model", par["model"] ]
    if par["total_droplets_included"]:
      cmd_pars += [ "--total-droplets-included", str(par["total_droplets_included"]) ]
    if par["epochs"]:
      cmd_pars += [ "--epochs", str(par["epochs"]) ]
    if par["fpr"]:
      cmd_pars += [ "--fpr", str(par["fpr"]) ]
    if par["exclude_antibody_capture"]:
      cmd_pars += [ "--exclude-antibody-capture" ]
    if par["learning_rate"]:
      cmd_pars += [ "--learning-rate", str(par["learning_rate"]) ]
    if par["cuda"]:
      cmd_pars += [ "--cuda" ]

    out = subprocess.check_output(cmd_pars).decode("utf-8")
    
    # have to use custom read_10x_h5 function for now
    # will be fixed when https://github.com/scverse/scanpy/pull/2344 is merged
    # adata_out = sc.read_10x_h5(output_file, gex_only=False)
    adata_out = anndata_from_h5(output_file)

    # store output data in mudata
    data.layers[par["layer_output"]] = adata_out.X

    obs_store = { 
      "obs_latent_rt_efficiency": "latent_RT_efficiency", 
      "obs_latent_cell_probability": "latent_cell_probability", 
      "obs_latent_scale": "latent_scale"
    }
    var_store = { 
      "var_ambient_expression": "ambient_expression"
    }
    obsm_store = { 
      "obsm_latent_gene_encoding": "latent_gene_encoding"
    }
    for to_name, from_name in obs_store.items():
      if par[to_name]:
        data.obs[par[to_name]] = adata_out.obs[from_name]
    for to_name, from_name in var_store.items():
      if par[to_name]:
        data.var[par[to_name]] = adata_out.var[from_name]
    for to_name, from_name in obsm_store.items():
      if par[to_name]:
        data.obsm[par[to_name]] = adata_out.obsm[from_name]


logger.info("Writing to file %s", par["output"])
mdata.write(filename=par["output"])