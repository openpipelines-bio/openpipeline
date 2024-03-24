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

# download foundational model files (full_human)
# https://drive.google.com/drive/folders/1oWh_-ZRdhtoGQ2Fw24HP41FgLoomVo-y
gdown '1H3E_MJ-Dl36AQV6jLbna2EdvgPaqvqcC' -O "${foundation_model_dir}/vocab.json"
gdown '1hh2zGKyWAx3DyovD30GStZ3QlzmSqdk1' -O "${foundation_model_dir}/args.json"
gdown '14AebJfGOUF047Eg40hk57HCtrb0fyDTm' -O "${foundation_model_dir}/best_model.pt"

# create test data dir
test_resources_dir="$OUT/test_resources"
mkdir -p "$test_resources_dir"

# download test data
# https://drive.google.com/file/d/1z_0vWYMhRuRiD1EyhuFtY9ReIR0msWaL/view?usp=drive_link
gdown '1z_0vWYMhRuRiD1EyhuFtY9ReIR0msWaL' -O "${test_resources_dir}/Kim2020_Lung.h5ad"

python <<HEREDOC
import anndata as ad
import mudata as mu
input_adata = ad.read_h5ad("${test_resources_dir}/Kim2020_Lung.h5ad")
input_mdata = mu.MuData({'rna': input_adata})
input_mdata.write_h5mu("${test_resources_dir}/Kim2020_Lung.h5mu")
HEREDOC

rm "${test_resources_dir}/Kim2020_Lung.h5ad"