#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparser"))
suppressPackageStartupMessages(library("data.table"))
parser <- arg_parser("Parameters for demuxlet summarzing")
parser <- add_argument(parser, "--demuxlet_out", help = "Path of demuxlet output file", default = NULL)
args <- parse_args(parser)
res <- fread(args$demuxlet_out)
res <- res[,c(2,5,6)]
res$`donor_id` <- sapply(res$BEST.GUESS,function(x){
    splitlist = strsplit(x,",")[[1]]
    if (splitlist[[1]] == splitlist[[2]]){
      splitlist[[2]]}
    else{
      "NA"}
})
res[`donor_id` =="NA",]$`donor_id` = "DBL"
res[DROPLET.TYPE=="AMB",]$`donor_id` = "AMB"
colnames(res)[1] <- "cell"
demuxlet_assign <- res[, c("cell", "donor_id")]

write.table(demuxlet_assign, paste0(dirname(args$demuxlet_out),"/assignment.tsv"), quote = FALSE)
