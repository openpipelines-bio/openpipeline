#!/bin/bash

set -eo pipefail

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

ID=scvi_tools
OUT=resources_test/$ID/$ID
DIR=$(dirname "$OUT")

# ideally, this would be a versioned pipeline run
[ -d "$DIR" ] || mkdir -p "$DIR"

# Prepare test data as in scvi-tools tutorial: https://docs.scvi-tools.org/en/stable/tutorials/notebooks/scarches_scvi_tools.html
python <<HEREDOC
import numpy as np
import scanpy as sc

# Download file with pancreas data
url = "https://figshare.com/ndownloader/files/24539828"
adata = sc.read("pancreas.h5ad", backup_url=url)

# Split file to reference and query by the technology used
query = np.array([s in ["smartseq2", "celseq2"] for s in adata.obs.tech])
adata_ref = adata[~query].copy()
adata_query = adata[query].copy()

# Subset HVG
sc.pp.highly_variable_genes(
    adata_ref,
    n_top_genes=2000,
    batch_key="tech",
    subset=True
)
adata_query = adata_query[:, adata_ref.var_names].copy()

sc.write("${OUT}_pancreas_query_test.h5ad", adata_query)
sc.write("${OUT}_pancreas_ref_test.h5ad", adata_ref)
HEREDOC

# convert 10x h5 to h5mu
bin/viash run src/convert/from_h5ad_to_h5mu/config.vsh.yaml -- \
  --input "${OUT}_pancreas_query_test.h5ad" \
  --output "${OUT}_pancreas_query_test.h5mu"

  bin/viash run src/convert/from_h5ad_to_h5mu/config.vsh.yaml -- \
  --input "${OUT}_pancreas_ref_test.h5ad" \
  --output "${OUT}_pancreas_ref_test.h5mu"

rm ${OUT}_pancreas_ref_test.h5ad
