requireNamespace("processx", quietly = TRUE)
requireNamespace("readr", quietly = TRUE)
library(dplyr, warn.conflicts = FALSE)

## VIASH START
par <- list(
  sam = "resources_test/demuxafy_test_data/pooled.sorted.bam",
  vcf = "resources_test/demuxafy_test_data/test_dataset.vcf",
  output = "freemuxlet_result",
  out = "out",
  nsample = "14"
)
## VIASH END

if (!dir.exists(par$output)) {
  dir.create(par$output, recursive = TRUE, showWarnings = FALSE)
}

cmd <- c(
  "popscle", "freemuxlet",
  "--out", paste0(par$output, "/", par$out)
)

argmap <- c(
  "plp" = "--plp",
  "init_cluster" = "--init-cluster",
  "nsample" = "--nsample",
  "verbose" = "--verbose",
  "doublet_prior" = "--doublet-prior",
  "geno_error" = "--geno-error",
  "bf_thres" = "--bf-thres",
  "frac_init_clust" = "--frac-init-clust",
  "iter_init" = "--iter-init",
  "seed" = "--seed",
  "cap_bq" = "--cap-BQ",
  "min_bq" = "--min-BQ",
  "min_total" = "--min-total",
  "min_umi" = "--min-umi",
  "min_snp" = "--min-snp",
  "group_list" = "--group-list",
  "aux_files" = "--aux-files",
  "keep_init_missing" = "--keep-init-missing",
  "randomize_singlet_score" = "randomize-singlet-score"
)

for (arg in names(argmap)) {
  if (!is.null(par[[arg]])) {
    if (arg %in% c("aux_files", "keep_init_missing", "randomize_singlet_score")) {
      if (toupper(par[[arg]]) == TRUE)
        cmd <- c(cmd, argmap[[arg]])
    }else {
      cmd <- c(cmd, argmap[[arg]], par[[arg]])
    }
  }
}

zzz <- processx::run(
  cmd[[1]],
  args = cmd[-1],
  echo = TRUE,
  echo_cmd = TRUE
)

if (zzz$status != 0) {
  stop("Command failed with status ", zzz$status)
}

out_file <- paste0(par$output, "/", par$out, ".clust1.samples.gz")
if (!file.exists(out_file)) {
  stop("Output file '", out_file, "' not found")
}

res <- readr::read_tsv(out_file)

res2 <- res %>%
  mutate(
    donor_part1 = gsub("([^,]*),([^,]*)*", "\\1", BEST.GUESS),
    donor_part2 = gsub("([^,]*),([^,]*)*", "\\2", BEST.GUESS),
    donor_id = case_when(
      donor_part1 == donor_part2 ~ donor_part1,
      TRUE ~ DROPLET.TYPE
    )
  )

freemuxlet_assign <- res2 %>% select(cell = BARCODE, donor_id)

readr::write_csv(
  freemuxlet_assign,
  paste0(par$output, "/cell_annotation.csv")
)