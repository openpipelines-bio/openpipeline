## VIASH START
par <- list(
  input = paste0(
    "resources_test/pbmc_1k_protein_v3/",
    "pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5"
  ),
  output = "output.h5",
  min_library_size = 1000,
  min_cells_per_gene = 300,
  keep_feature_types = "Antibody Capture",
  verbose = TRUE
)
## VIASH END

if (par$verbose) cat("Loading dependencies\n")
requireNamespace("hdf5r", quietly = TRUE)

if (par$verbose) cat("Opening h5 file\n")
h5 <- hdf5r::H5File$new(par$input, mode = "r")

if (par$verbose) cat("Reading data in memory\n")
features__all_tag_keys <- h5[["matrix/features/_all_tag_keys"]][]

features <- data.frame(
  feature_type = h5[["matrix/features/feature_type"]][],
  genome = h5[["matrix/features/genome"]][],
  id = h5[["matrix/features/id"]][],
  name = h5[["matrix/features/name"]][]
)

mat <- Matrix::sparseMatrix(
  i = h5[["matrix/indices"]][],
  p = h5[["matrix/indptr"]][],
  x = h5[["matrix/data"]][],
  dims = h5[["matrix/shape"]][],
  index1 = FALSE,
  dimnames = list(
    features$id,
    h5[["matrix/barcodes"]][]
  )
)

if (par$verbose) {
  cat("Filtering out cells with library size < ",
    par$min_library_size, "\n",
    sep = ""
  )
}
library_size <- Matrix::colSums(mat)
mat2 <- mat[, library_size >= par$min_library_size, drop = FALSE]

if (par$verbose) {
  cat("Filtering genes with num cells < ",
    par$min_cells_per_gene, "\n",
    sep = ""
  )
}
num_cells <- Matrix::rowSums(mat2 > 0)
mat3 <- mat2[
  num_cells >= par$min_cells_per_gene |
    features$feature_type %in% par$keep_feature_types, ,
  drop = FALSE
]
features2 <- features[match(rownames(mat3), features$id), , drop = FALSE]

# helper fun
set_with_type <- function(path, value) {
  orig_dtype <- h5[[path]]$get_type()
  orig_chunk <- h5[[path]]$chunk_dims
  if (is.na(orig_chunk)) orig_chunk <- "auto"
  h5new$create_dataset(path, value, dtype = orig_dtype, chunk_dims = orig_chunk)
}

# create new file
if (par$verbose) cat("Saving h5 file at '", par$output, "'\n", sep = "")
h5new <- hdf5r::H5File$new(par$output, mode = "w")
zz <- h5new$create_group("matrix")
zz <- h5new$create_group("matrix/features")

set_with_type("matrix/features/feature_type", features2$feature_type)
set_with_type("matrix/features/genome", features2$genome)
set_with_type("matrix/features/id", features2$id)
set_with_type("matrix/features/name", features2$name)
set_with_type("matrix/features/_all_tag_keys", features__all_tag_keys)
set_with_type("matrix/indices", mat3@i)
set_with_type("matrix/indptr", mat3@p)
set_with_type("matrix/data", as.integer(mat3@x))
set_with_type("matrix/shape", dim(mat3))
set_with_type("matrix/barcodes", colnames(mat3))

for (attname in hdf5r::h5attr_names(h5)) {
  h5new$create_attr(attname, hdf5r::h5attr(h5, attname))
}
h5new$close_all()
h5$close_all()
