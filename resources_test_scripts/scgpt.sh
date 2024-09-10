set -eo pipefail

# ensure that the command below is run from the root of the repository
REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

# settings
ID=scgpt
OUT=resources_test/$ID

# create foundational model directory
foundation_model_dir="$OUT/source"
mkdir -p "$foundation_model_dir"

# install gdown if necessary
# Check whether gdown is available
if ! command -v gdown &> /dev/null; then
    echo "This script requires gdown. Please make sure the binary is added to your PATH."
    exit 1
fi

echo "> Downloading scGPT foundation model (full_human)"
# download foundational model files (full_human)
# https://drive.google.com/drive/folders/1oWh_-ZRdhtoGQ2Fw24HP41FgLoomVo-y
gdown '1H3E_MJ-Dl36AQV6jLbna2EdvgPaqvqcC' -O "${foundation_model_dir}/vocab.json"
gdown '1hh2zGKyWAx3DyovD30GStZ3QlzmSqdk1' -O "${foundation_model_dir}/args.json"
gdown '14AebJfGOUF047Eg40hk57HCtrb0fyDTm' -O "${foundation_model_dir}/best_model.pt"

# create test data dir
test_resources_dir="$OUT/test_resources"
mkdir -p "$test_resources_dir"

echo "> Downloading test resources"
# download test data
# https://drive.google.com/file/d/1z_0vWYMhRuRiD1EyhuFtY9ReIR0msWaL/view?usp=drive_link
gdown '1z_0vWYMhRuRiD1EyhuFtY9ReIR0msWaL' -O "${test_resources_dir}/Kim2020_Lung.h5ad"

echo "> Converting to h5mu"
python <<HEREDOC
import anndata as ad
import mudata as mu
input_adata = ad.read_h5ad("${test_resources_dir}/Kim2020_Lung.h5ad")
input_mdata = mu.MuData({'rna': input_adata})
input_mdata.write_h5mu("${test_resources_dir}/Kim2020_Lung.h5mu")
HEREDOC

echo "> Subsetting datasets"
viash run src/filter/subset_h5mu/config.vsh.yaml -p docker -- \
  --input "${test_resources_dir}/Kim2020_Lung.h5mu" \
  --output "${test_resources_dir}/Kim2020_Lung_subset.h5mu" \
  --number_of_observations 4000

rm "${test_resources_dir}/Kim2020_Lung.h5ad"

echo "> Preprocessing datasets"
nextflow \
  run . \
  -main-script target/nextflow/workflows/multiomics/process_samples/main.nf \
  -profile docker \
  -c src/workflows/utils/labels_ci.config \
  --input "${test_resources_dir}/Kim2020_Lung_subset.h5mu" \
  --output "Kim2020_Lung_subset_preprocessed.h5mu" \
  --publish_dir "${test_resources_dir}"

echo "> Filtering highly variable features"
viash run src/feature_annotation/highly_variable_features_scanpy/config.vsh.yaml -p docker -- \
  --input "${test_resources_dir}/iKim2020_Lung_subset_preprocessed.h5mu" \
  --output "${test_resources_dir}/Kim2020_Lung_subset_hvg.h5mu" \
  --layer "log_normalized" \
  --var_name_filter "filter_with_hvg" \
  --n_top_features 1200 \
  --flavor "seurat_v3"

viash run src/filter/do_filter/config.vsh.yaml -p docker -- \
  --input "${test_resources_dir}/Kim2020_Lung_subset_hvg.h5mu" \
  --output "${test_resources_dir}/Kim2020_Lung_subset_hvg_filtered.h5mu" \
  --var_filter "filter_with_hvg"
  
echo "> Running scGPT cross check genes"
viash run src/scgpt/cross_check_genes/config.vsh.yaml -p docker -- \
  --input "${test_resources_dir}/Kim2020_Lung_subset_hvg_filtered.h5mu" \
  --output "${test_resources_dir}/Kim2020_Lung_subset_genes_cross_checked.h5mu" \
  --vocab_file "${foundation_model_dir}/vocab.json"

echo "> Running scGPT binning"
viash run src/scgpt/binning/config.vsh.yaml -p docker -- \
  --input "${test_resources_dir}/Kim2020_Lung_subset_genes_cross_checked.h5mu" \
  --input_layer "log_normalized" \
  --output "${test_resources_dir}/Kim2020_Lung_subset_binned.h5mu"

echo "> Running scGPT tokenizing"
viash run src/scgpt/pad_tokenize/config.vsh.yaml -p docker -- \
  --input "${test_resources_dir}/Kim2020_Lung_subset_binned.h5mu" \
  --input_layer "binned" \
  --output "${test_resources_dir}/Kim2020_Lung_subset_tokenized.h5mu" \
  --model_vocab "${foundation_model_dir}/vocab.json"

echo "> Running scGPT integration"
viash run src/scgpt/embedding/config.vsh.yaml -p docker -- \
  --input "${test_resources_dir}/Kim2020_Lung_subset_tokenized.h5mu" \
  --output "${test_resources_dir}/Kim2020_Lung_subset_scgpt_integrated.h5mu" \
  --model "${foundation_model_dir}/best_model.pt" \
  --model_vocab "${foundation_model_dir}/vocab.json" \
  --model_config "${foundation_model_dir}/args.json" \
  --obs_batch_label "sample"

echo "> Removing unnecessary files in test resources dir"
find "${test_resources_dir}" -type f \( ! -name "Kim2020_*" -o ! -name "*.h5mu" \) -delete

echo "> scGPT test resources are ready!"

