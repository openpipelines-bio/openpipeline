#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparser"))
suppressPackageStartupMessages(library("data.table"))
parser <- arg_parser("Parameters for freemuxlet summarzing")
parser <- add_argument(parser, "--freemuxlet_out",
            help = "Path of freemuxlet output file", default = NULL)
args <- parse_args(parser)
res <- fread(args$freemuxlet_out)
res <- res[, c(2, 5, 6)]
res$`donor_id` <- sapply(res$BEST.GUESS, function(x) {
    splitlist <- strsplit(x,",")[[1]]
    if (splitlist[[1]] == splitlist[[2]]) {
      splitlist[[2]]
      }
      else{
        "NA"
      }
})
res[`donor_id` == "NA", ]$`donor_id` <- "doublet"
res[DROPLET.TYPE == "AMB", ]$`donor_id` <- "unassigned"
colnames(res)[1] <- "cell"
freemuxlet_assign <- res[, c("cell", "donor_id")]

write.csv(freemuxlet_assign,
    paste0(dirname(args$freemuxlet_out), "/cell_annotation.csv"),
    row.names = FALSE, quote = FALSE)
