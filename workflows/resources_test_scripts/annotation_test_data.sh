#!/bin/bash

set -eo pipefail

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

ID=annotation_test_data
OUT=resources_test/$ID/
OUT_ONTOLOGY=resources_test/$ID/popv_cl_ontology/
DIR=$(dirname "$OUT")
DIR_ONTOLOGY=$(dirname "$OUT_ONTOLOGY")

# ideally, this would be a versioned pipeline run
[ -d "$DIR" ] || mkdir -p "$DIR"
[ -d "$DIR_ONTOLOGY" ] || mkdir -p "$DIR_ONTOLOGY"

# dataset page:
# https://www.10xgenomics.com/resources/datasets/1-k-pbm-cs-from-a-healthy-donor-gene-expression-and-cell-surface-protein-3-standard-3-0-0


# Download query h5
wget https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5 \
  -O "${OUT}/tmp_blood_test_query.h5"


# # Download Tabula Sapiens Blood reference h5ad
wget https://www.dropbox.com/s/4cg6zj340oelhlg/Blood.h5ad?dl=1 \
  -O "${OUT}/tmp_blood_test_reference.h5ad"

# Download PopV specific CL ontology files
wget https://raw.githubusercontent.com/czbiohub/PopV/main/ontology/cl.obo \
-O "${OUT_ONTOLOGY}/cl.obo"
wget https://raw.githubusercontent.com/czbiohub/PopV/main/ontology/cl.ontology \
-O "${OUT_ONTOLOGY}/cl.ontology"
wget https://raw.githubusercontent.com/czbiohub/PopV/main/ontology/cl.ontology.nlp.emb \
-O "${OUT_ONTOLOGY}/cl.ontology.nlp.emb"


# Local python MUST have
#   scanpy==1.8.2
#   anndata==0.7.8
python <<HEREDOC
import scanpy as sc
import anndata
assert sc.__version__ == '1.8.2', 'Please make sure to have scanpy 1.8.2 version installed (due to requirement PopV)!'
assert anndata.__version__ == '0.7.8', 'Please make sure to have anndata 0.7.8 version installed (due to requirement PopV)!'
HEREDOC


# Process query h5 
# Convert h5 to h5ad and comply format
python <<HEREDOC
import scanpy as sc
import anndata
query_adata = sc.read_10x_h5("${OUT}/tmp_blood_test_query.h5", gex_only=False)
query_adata.var_names_make_unique()
gex_query_adata = query_adata[:, query_adata.var["feature_types"] == "Gene Expression"]
gex_query_adata.obs['batch'] = 'blood_test_query_sample'
gex_query_adata.raw = gex_query_adata
assert gex_query_adata.shape == (713, 33538)
gex_query_adata.write("${OUT}/blood_test_query.h5ad", compression='gzip')
HEREDOC
rm "${OUT}/tmp_blood_test_query.h5"


# Process Tabula Sapiens Blood reference h5ad
# Select one individual and 100 cells per cell type)
python <<HEREDOC
import scanpy as sc
import anndata
ref_adata = sc.read_h5ad("${OUT}/tmp_blood_test_reference.h5ad")
sub_ref_adata = ref_adata[ref_adata.obs["donor_method"] == "TSP1410X"] 
n=100
s=sub_ref_adata.obs.groupby('cell_ontology_class').cell_ontology_class.transform('count')
sub_ref_adata_final = sub_ref_adata[sub_ref_adata.obs[s>=n].groupby('cell_ontology_class').head(n).index]
assert sub_ref_adata_final.shape == (500, 58870)
sub_ref_adata_final.write("${OUT}/blood_test_reference.h5ad", compression='gzip')
HEREDOC
rm "${OUT}/tmp_blood_test_reference.h5ad"
