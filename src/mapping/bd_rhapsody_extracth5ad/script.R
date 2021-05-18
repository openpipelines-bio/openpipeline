## VIASH START
par <- list(
  input = "/home/rcannood/workspace/jnj/bdrhap/out2/",
  output = "/home/rcannood/workspace/jnj/bdrhap/out2/output.h5ad"
)
## VIASH END

library(tidyverse)
library(anndata)

# detect file
metrics_summary <- list.files(par$input, pattern = "Metrics_Summary.csv$", full.names = TRUE)
reads_per_cell <- list.files(par$input, pattern = "ReadsPerCell_Unfiltered.csv.gz$", full.names = TRUE) %>% readr::read_csv(skip = 6)

rpc <- reads_per_cell %>%
  tibble::column_to_rownames("Cell_Index") %>%
  as.matrix %>%
  Matrix::Matrix(sparse = TRUE)

ad <- anndata::AnnData(
  X = rpc
)

ad$raw <- ad
ad$var["feature_types"] <- "Gene Expression"

# todo: can we add more metadata?

ad$write_h5ad(filename = par$output, compression = "gzip")
