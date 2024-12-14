library(dplyr)
library(purrr)

## VIASH START
par <- list(
  base_dir = "src",
  pattern = "*.vsh.yaml",
  n_dirname_drop = 1,
  n_basename_id = 1,
  output = "output.yaml",
  id_name = "id",
  path_name = "path",
  group_name = "param_list"
)
## VIASH END

cat("> Listing files of base dir ", par$base_dir, "\n", sep = "")
paths <- list.files(
  normalizePath(par$base_dir),
  pattern = par$pattern,
  recursive = TRUE,
  full.names = TRUE
)

cat("> Traversing up ", par$n_dirname_apply, " times\n", sep = "")
for (i in seq_len(par$n_dirname_drop)) {
  paths <- dirname(paths) %>% unique()
}

# removing /viash_automount in case we're inside a docker container
paths <- gsub("^/viash_automount", "", paths)

cat("> Checking whether basenames are unique\n")
i <- par$n_basename_id
maxi <- strsplit(paths, "/") %>%
  map_int(length) %>%
  max()

regex <- paste0(".*/(", paste(rep("[^/]+/", i), collapse = ""), "[^/]*)$")
ids <- gsub("/", "_", gsub(regex, "\\1", paths))

cat("> Printing first five rows\n")
print(tibble(id = ids, path = paths) %>% head(5))
cat("\n")

while (i < maxi && any(duplicated(ids))) {
  i <- i + 1
  cat(
    "Duplicated ids detected, combining with ", i,
    " dirnames in an attempt to get unique ids.\n"
  )
  regex <- paste0(".*/(", paste(rep("[^/]+/", i), collapse = ""), "[^/]*)$")
  ids <- gsub("/", "_", gsub(regex, "\\1", paths))

  cat("> Printing first five rows\n")
  print(tibble(id = ids, path = paths) %>% head(5))
  cat("\n")
}

cat("> Transforming into list of items\n")
par_list <- map2(
  ids, paths,
  function(id, input) {
    setNames(list(id, input), c(par$id_name, par$path_name))
  }
)

if (!is.null(par$group_name)) {
  par_list <- setNames(list(par_list), par$group_name)
}

cat("> Writing as YAML\n")
yaml::write_yaml(par_list, par$output)
