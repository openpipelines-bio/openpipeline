#!/bin/bash

set -eo pipefail

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

ID=annotation_test_data
OUT=resources_test/$ID/
DIR=$(dirname "$OUT")

# ideally, this would be a versioned pipeline run
[ -d "$DIR" ] || mkdir -p "$DIR"

# Download Tabula Sapiens Blood reference h5ad from https://doi.org/10.5281/zenodo.7587774
# wget "https://zenodo.org/record/7587774/files/TS_Blood_filtered.h5ad?download=1" -O "${OUT}/tmp_TS_Blood_filtered.h5ad"

# # Download Tabula Sapiens Blood pretrained model from https://doi.org/10.5281/zenodo.7580707
# wget "https://zenodo.org/record/7580707/files/pretrained_models_Blood_ts.tar.gz?download=1" -O "${OUT}/tmp_pretrained_models_Blood_ts.tar.gz"

# # # Download PopV specific CL ontology files - now needed for OnClass
# # OUT_ONTOLOGY="${OUT}/ontology"
# # [ -d "$OUT_ONTOLOGY" ] || mkdir -p "$OUT_ONTOLOGY"
# # wget https://raw.githubusercontent.com/czbiohub/PopV/main/ontology/cl.obo \
# # -O "${OUT_ONTOLOGY}/cl.obo"
# # wget https://raw.githubusercontent.com/czbiohub/PopV/main/ontology/cl.ontology \
# # -O "${OUT_ONTOLOGY}/cl.ontology"
# # wget https://raw.githubusercontent.com/czbiohub/PopV/main/ontology/cl.ontology.nlp.emb \
# # -O "${OUT_ONTOLOGY}/cl.ontology.nlp.emb"


# # Process Tabula Sapiens Blood reference h5ad
# # (Select one individual and 100 cells per cell type)
# python <<HEREDOC
# import anndata as ad
# ref_adata = ad.read_h5ad("${OUT}/tmp_TS_Blood_filtered.h5ad")
# sub_ref_adata = ref_adata[ref_adata.obs["donor_assay"] == "TSP14_10x 3' v3"] 
# n=100
# s=sub_ref_adata.obs.groupby('cell_ontology_class').cell_ontology_class.transform('count')
# sub_ref_adata_final = sub_ref_adata[sub_ref_adata.obs[s>=n].groupby('cell_ontology_class').head(n).index]
# # assert sub_ref_adata_final.shape == (500, 58870)
# sub_ref_adata_final.write("${OUT}/TS_Blood_filtered.h5ad", compression='gzip')
# HEREDOC

# # rm "${OUT}/tmp_blood_test_reference.h5ad"

# echo "Fetching OnClass data and models"

# OUT_ONTOLOGY="${OUT}/ontology"
# [ -d "$OUT_ONTOLOGY" ] || mkdir -p "$OUT_ONTOLOGY"
# wget https://figshare.com/ndownloader/files/28394466 -O "${OUT_ONTOLOGY}/OnClass_data_public_minimal.tar.gz"
# tar -xzvf "${OUT_ONTOLOGY}/OnClass_data_public_minimal.tar.gz" -C "${OUT_ONTOLOGY}" --strip-components=2
# rm "${OUT_ONTOLOGY}/allen.ontology"
# rm "${OUT_ONTOLOGY}/OnClass_data_public_minimal.tar.gz"

wget https://figshare.com/ndownloader/files/28394541 -O "${OUT}/OnClass_models.tar.gz"
tar -xzvf "${OUT}/OnClass_models.tar.gz" -C "${OUT}" --strip-components=1
rm "${OUT}/OnClass_models.tar.gz"

find "${OUT}/Pretrained_model" ! -name "example_file_model*" -type f -exec rm -f {} +
mv "${OUT}/Pretrained_model" "${OUT}/onclass_model"
