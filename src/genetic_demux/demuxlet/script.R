requireNamespace("processx", quietly = TRUE)
requireNamespace("readr", quietly = TRUE)
library(dplyr, warn.conflicts = FALSE)

## VIASH START
par <- list(
  sam = "resources_test/demuxafy_test_data/chr_1_pooled.sorted.bam",
  vcf = "resources_test/demuxafy_test_data/test_dataset.vcf",
  output = "demuxlet_result",
  out = "out",
  field = "GP"
)
## VIASH END

if (!dir.exists(par$output)) {
  dir.create(par$output, recursive = TRUE, showWarnings = FALSE)
}

cmd <- c(
  "popscle", "demuxlet",
  "--out", paste0(par$output, "/", par$out)
)

argmap <- c(
  "tag_group" = "--tag-group",
  "tag_umi" = "--tag-UMI",
  "field" = "--field",
  "geno_error_offset" = "--geno-error-offset",
  "geno_error_coeff" = "--geno-error-coeff",
  "r2_info" = "--r2-info",
  "min_mac" = "--min-mac",
  "min_call_rate" = "--min-callrate",
  "alpha" = "--alpha",
  "doublet_prior" = "--doublet-prior",
  "sm" = "--sm",
  "sm_list" = "--sm-list",
  "sam_verbose" = "--sam-verbose",
  "vcf_verbose" = "--vcf-verbose",
  "cap_bq" = "--cap-BQ",
  "min_bq" = "--min-BQ",
  "min_mq" = "--min-MQ",
  "min_td" = "--min-TD",
  "excl_flag" = "--excl-flag",
  "group_list" = "--group-list",
  "min_total" = "--min-total",
  "min_snp" = "--min-snp",
  "min_umi" = "--min-umi",
  "plp" = "--plp",
  "vcf" = "--vcf",
  "sam" = "--sam",
  "sm" = "--sm",
  "sm_list" = "--sm-list",
  "group_list" = "--group-list"
)

for (arg in names(argmap)) {
  if (!is.null(par[[arg]])) {
    cmd <- c(cmd, argmap[[arg]], par[[arg]])
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

out_file <- paste0(par$output, "/", par$out, ".best")
if (!file.exists(out_file)) {
  stop("Output file '", out_file, "' not found")
}
res <- readr::read_tsv(out_file)

res2 <- res %>%
  mutate(
    donor_part1 = gsub("([^,]*),([^,]*),.*", "\\1", BEST.GUESS),
    donor_part2 = gsub("([^,]*),([^,]*),.*", "\\2", BEST.GUESS),
    donor_id = case_when(
      donor_part1 == donor_part2 ~ donor_part1,
      TRUE ~ DROPLET.TYPE
    )
  )

demuxlet_assign <- res2 %>% select(cell = BARCODE, donor_id)

readr::write_csv(
  demuxlet_assign,
  paste0(par$output, "/cell_annotation.csv")
)
