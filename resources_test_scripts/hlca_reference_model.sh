#!/bin/bash

set -eo pipefail

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

ID=HLCA_reference_model
OUT=resources_test/$ID/$ID
DIR=$(dirname "$OUT")

# ideally, this would be a versioned pipeline run
[ -d "$DIR" ] || mkdir -p "$DIR"

# download and unarchive pre-trained scANVI model
wget https://zenodo.org/record/6337966/files/HLCA_reference_model.zip \
  -O "${OUT}.zip"

# # Test query data
# # Source publication: Delorey, Toni M., et al. “COVID-19 tissue atlases reveal SARS-CoV-2 pathology and cellular targets.” Nature 595.7865 (2021): 107-113.
# wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5230nnn/GSM5230027/suppl/GSM5230027_04-P103142-S149-R01_raw_feature_bc_matrix.h5.gz \
#   -O "${OUT}_query_test.h5.gz"
# gzip -d "${OUT}_query_test.h5.gz"

# # Prepare test data as in scvi-tools tutorial: https://docs.scvi-tools.org/en/stable/tutorials/notebooks/query_hlca_knn.html
# python <<HEREDOC
# import pandas as pd
# import scanpy as sc

# geo_metadata_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE171nnn/GSE171668/suppl/GSE171668_lung_metadata.csv.gz"
# metadata = pd.read_csv(geo_metadata_url, index_col=0)

# DATA_PATH = "${OUT}_query_test.h5"
# query_data = sc.read_10x_h5(DATA_PATH)
# # clean up .var.index (gene names)
# query_data.var['gene_names'] = query_data.var.index
# query_data.var.index = [idx.split("___")[-1] for idx in query_data.var.gene_ids]
# # clean up cell barcodes:
# query_data.obs.index = query_data.obs.index.str.rstrip("-1")
# # read in metadata (to select only cells of interest and remove empty drops)
# # subset to cells from our sample
# metadata = metadata.loc[metadata.donor == "D12_4",:].copy()
# # clean up barcodes:
# metadata.index = [idx.split("-")[-1] for idx in metadata.index]
# # subset adata to cells in metadata:
# query_data = query_data[metadata.index,:].copy()
# # add dataset information:
# query_data.obs['dataset'] = "test_dataset_delorey_regev"
# sc.write(DATA_PATH, query_data)
# HEREDOC

# # convert 10x h5 to h5mu
# viash run src/convert/from_h5ad_to_h5mu/config.vsh.yaml -- \
#   --input "${OUT}_query_test.h5" \
#   --output "${OUT}_query_test.h5mu"
