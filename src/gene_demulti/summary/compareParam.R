#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("R.utils"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("ggpubr"))


# transform analysis csv --------------------------------------------------------------------------------------------------------------------------------------------------------
transform_df <- function(result_melt){
  colnames(result_melt) = c('Trial','donor_identity','number')
  result_melt_sorted <- arrange(result_melt, Trial, number)
  result_melt_sorted_cumsum <- ddply(result_melt_sorted, "Trial", transform, label_ypos=cumsum(number))
  result_melt_sorted_cumsum <- result_melt_sorted_cumsum %>% 
  group_by('donor_identity') %>%
  dplyr::mutate(label_ypos_new = mean(label_ypos)) %>% 
  ungroup() %>% 
  as.data.frame()
  return (result_melt_sorted_cumsum)
}


# plot analysis csv ----------------------------------------------------------------------------------------------------------------------------------------------------------
bar_plot_group <- function(result_melt, tool){
  max_ylim = max(result_melt$label_ypos) + 50
  
  if(nrow(result_melt[result_melt$number==0,]) != 0){
    result_melt[result_melt$number==0,]$label_ypos_new = NA
    result_melt[result_melt$number==0,]$label_ypos = NA
  }
  
  bar_plot_g = ggplot(result_melt, aes(x = Trial, y = number, fill = donor_identity)) + 
    geom_bar(stat="identity") + 
    geom_text(aes(y = label_ypos, label = number), vjust = 2, size = 3.5) +
    scale_fill_brewer(palette="Blues") +ylim(0, max_ylim) + ggtitle(paste0("Demultiplexing Results: ", tool)) + 
    theme_light() + theme(plot.title = element_text(hjust = 0.5))
    
  png(paste0(tool,"_barplotgroup.png"))
  print(bar_plot_g)
  while (!is.null(dev.list()))  dev.off() 
}


# read demuxlet file ----------------------------------------------------------------------------------------------------------------------------------------------------------
readDemux <- function(demux){
  demux <- fread(demux)
  demux_df <- demux[,c(2,5,6,7)]
  demux_df$donor.identity  <- sapply(demux_df$BEST.GUESS,function(x){
    splitlist = strsplit(x,",")[[1]]
    if (splitlist[[1]] == splitlist[[2]]){ 
      splitlist[[2]]}
    else{
      "NA"}
  })
  demux_df[donor.identity =="NA",]$donor.identity = "DBL"
  demux_df[DROPLET.TYPE=="AMB",]$donor.identity = "AMB"
  return (demux_df)
}

# analysis demuxlet data ----------------------------------------------------------------------------------------------------------------------------------------------------------

num_stat_demux <- function(demux_list, barcode){
  result_df = as.data.frame(matrix(nrow= length(demux_list), ncol=4))
  colnames(result_df) = c("Trial","SNG","DBL","AMB")
  
  f1 = readDemux(demux_list[1])
  
  donor_identity = sort(unique(append(f1$donor.identity, c("AMB","DBL"))))
  stat_donor_identity = as.data.frame(matrix(nrow= length(demux_list), ncol=length(donor_identity)+1))
  colnames(stat_donor_identity) = append("Trial",donor_identity)
  
  num_barcode = length(barcode)
  assignment = as.data.frame(matrix(ncol= 1, nrow=num_barcode))
  assignment[,1] = barcode
  colnames(assignment) = c("barcode")
  
  trial = as.data.frame(matrix(nrow=29, ncol=1))
  colnames(trial) = c("Parameter")
  trial[,1] = c("bam","tag-group","tag-UMI","vcf-ref", "sm","sm-list","sam-verbose","vcf-verbose","skip-umi", "cap-BQ","min-BQ","min-MQ","min-TD","excl-flag", "group-list","min-total", "min-uniq", "min-umi","min-snp", "plp", "vcf-donor", "field","geno-error-offset","geno-error-coeff","r2-info","min-mac","min-callrate","alpha","doublet-prior")
  
  rindex = 1

  for (f in demux_list){
    demux_df = readDemux(f)
    num_sin = nrow(demux_df[demux_df$DROPLET.TYPE=="SNG",])
    num_dou = nrow(demux_df[demux_df$DROPLET.TYPE=="DBL",])
    num_amb = nrow(demux_df[demux_df$DROPLET.TYPE=="AMB",])

    row_name = paste("Demuxlet Trial",rindex)
    result_df[rindex,1] = row_name
    result_df[rindex,2] = num_sin
    result_df[rindex,3] = num_dou
    result_df[rindex,4] = num_amb
    
    donor_identity = unique(demux_df$donor.identity)
    for(c in donor_identity){
      num_cell = nrow(demux_df[demux_df$donor.identity==c,])
      stat_donor_identity[rindex,c] = num_cell
    }
    stat_donor_identity[rindex,1] = row_name
    stat_donor_identity[is.na(stat_donor_identity)] = 0
    
    assignment[,row_name] = merge(assignment, demux_df, by.x = "barcode", by.y = "BARCODE", all = TRUE)$donor.identity
    
    f = strsplit(f, split='/demux+',fixed=TRUE)
    if(length(f[[1]]) == 1){
        f = f[[1]][1]
    }
    else{
	      f = f[[1]][2]
    }
    trial_param = strsplit(substr(f,0,nchar(f)-5), '+',fixed = TRUE)[[1]]
    trial[,paste("Demuxlet Trial",rindex)] = trial_param 
    
    rindex = rindex + 1
  }
  assignment = assignment[order(assignment$barcode),]
  write.table(result_df, paste0(outdir, "demuxlet_result.csv"), row.names=FALSE, quote = FALSE, sep = ';')
  write.table(stat_donor_identity, paste0(outdir, "demuxlet_result_detail.csv"), row.names=FALSE, quote = FALSE, sep = ';')
  write.table(assignment, paste0(outdir, "demuxlet_assignment.csv"), row.names=FALSE, quote = FALSE, sep = ';')
  write.table(trial, paste0(outdir, "demuxlet_trial.csv"), row.names=FALSE, quote = FALSE, sep = ';')
  
}
# analysis freemuxlet data -------------------------------------------------------------------------------------------------

num_stat_freemuxlet <- function(freemuxlet_list, barcode){
  result_df = as.data.frame(matrix(nrow= length(freemuxlet_list), ncol=4))
  colnames(result_df) = c("Trial","SNG","DBL","AMB")
  f1 = readDemux(freemuxlet_list[1])
  
  donor_identity = sort(unique(append(f1$donor.identity, c("AMB","DBL"))))
  stat_donor_identity = as.data.frame(matrix(nrow= length(freemuxlet_list), ncol=length(donor_identity)+1))
  colnames(stat_donor_identity) = append("Trial",donor_identity)
  
  num_barcode = length(barcode)
  assignment = as.data.frame(matrix(ncol= 1, nrow=num_barcode))
  assignment[,1] = barcode
  colnames(assignment) = c("barcode")
  
  trial = as.data.frame(matrix(nrow=29, ncol=1))
  colnames(trial) = c("Parameter")
  trial[,1] = c("bam","tag-group","tag-UMI","vcf-ref","sm","sm-list","sam-verbose","vcf-verbose","skip-umi","cap-BQ","min-BQ","min-MQ", "min-TD","excl-flag","group-list","min-total", "min-uniq","min-umi", "min-snp","plp","init-cluster", "nsample", "aux-files", "verbose", "doublet-prior", "bf-thres", "frac-init-clust", "iter-init", "keep-init-missing")
  
  rindex = 1
  
  for (f in freemuxlet_list){
    freemuxlet_df = readDemux(f)
    num_sin = nrow(freemuxlet_df[freemuxlet_df$DROPLET.TYPE=="SNG",])
    num_dou = nrow(freemuxlet_df[freemuxlet_df$DROPLET.TYPE=="DBL",])
    num_amb = nrow(freemuxlet_df[freemuxlet_df$DROPLET.TYPE=="AMB",])
    row_name = paste("Freemuxlet Trial",rindex)
    result_df[rindex,1] = row_name
    result_df[rindex,2] = num_sin
    result_df[rindex,3] = num_dou
    result_df[rindex,4] = num_amb
    
    donor_identity = unique(freemuxlet_df$donor.identity)
    for(c in donor_identity){
      num_cell = nrow(freemuxlet_df[freemuxlet_df$donor.identity==c,])
      stat_donor_identity[rindex,c] =num_cell
    }
    stat_donor_identity[rindex,1] = row_name
    stat_donor_identity[is.na(stat_donor_identity)] = 0
    
    assignment[,row_name] = merge(assignment, freemuxlet_df, by.x = "barcode", by.y = "BARCODE", all = TRUE)$donor.identity
  
    f = strsplit(f, split='/freemuxlet+',fixed=TRUE)
    if(length(f[[1]]) == 1){
      f = f[[1]][1]
    }else{
      f = f[[1]][2]
    }
    trial_param = strsplit(substr(f,0,nchar(f)-29), '+',fixed = TRUE)[[1]]
    trial[,paste("Freemuxlet Trial",rindex)] = trial_param 
    
    rindex = rindex + 1
  }
  assignment = assignment[order(assignment$barcode),]
  write.table(result_df, paste0(outdir, "freemuxlet_result.csv"), row.names=FALSE, quote = FALSE, sep = ';')
  write.table(stat_donor_identity, paste0(outdir, "freemuxlet_result_detail.csv"), row.names=FALSE, quote = FALSE, sep = ';')
  write.table(assignment, paste0(outdir, "freemuxlet_assignment.csv"), row.names=FALSE, quote = FALSE, sep = ';')
  write.table(trial, paste0(outdir, "freemuxlet_trial.csv"), row.names=FALSE, quote = FALSE, sep = ';')
  
}


# analysis vireo data ----------------------------------------------------------------------------------------------------------------------------------------------------------

num_stat_vireo <- function(vireo_list, barcode){
  result_df = as.data.frame(matrix(nrow= length(vireo_list), ncol=4))
  colnames(result_df) = c("Trial","SNG","DBL","AMB")

  stat_donors = as.data.frame(matrix(nrow = 1))
  colnames(stat_donors) = "Trial"
  
  num_barcode = length(barcode)
  assignment = as.data.frame(matrix(ncol= 1, nrow=num_barcode))
  assignment[,1] = barcode
  colnames(assignment) = c("barcode")
  
  doublet_info = as.data.frame(matrix(ncol= 1, nrow=num_barcode))
  doublet_info[,1] = barcode
  colnames(doublet_info) = c("barcode")
  
  trial = as.data.frame(matrix(nrow=16, ncol=1))
  colnames(trial) = c("Paramter")
  trial[,1] = c("celldata","n-donor","vartrix-data","donor-file","geno-tag","no-doublet","n-init","extra-donor",
                       "extra-donor-mode","force-learn-GT","ase-mode","no-plot","random-seed",
                       "cell-range","call-ambient-rna","n-proc")
  
  rindex = 1
  
  for (f in vireo_list){
    summary = fread(paste(f,"/summary.tsv",sep = ""))
    donorsid = fread(paste(f,"/donor_ids.tsv",sep = ""))
    
    summary[summary == "doublet"] = "DBL"
    summary[summary == "unassigned"] = "AMB"
    donorsid[donorsid == "doublet"] = "DBL"
    donorsid[donorsid == "unassigned"] = "AMB"
    
    num_dou = summary[summary$Var1 == 'DBL',]$Freq[1]
    num_amb = summary[summary$Var1 == 'AMB',]$Freq[1]
    num_drop = nrow(summary)
    num_sin = num_drop - num_dou - num_amb
    
    row_name = paste("Vireo Trial",rindex)
    result_df[rindex,1] = row_name
    result_df[rindex,2] = num_sin
    result_df[rindex,3] = num_dou
    result_df[rindex,4] = num_amb
    
    donorsid$donor_id = str_replace_all(donorsid$donor_id,"donor","")
    donorsid$best_doublet = str_replace_all(donorsid$best_doublet,"donor","")
    
    doublet_info[,c(row_name)] = NA
    doublet_list = donorsid[donorsid$donor_id %in% c("AMB","DBL"),c('cell','best_doublet')]
    doublet_info[doublet_info$barcode %in% doublet_list$cell,c(row_name)] = doublet_list[,"best_doublet"]

    summary$Var1 = str_replace_all(summary$Var1,"donor","")
    
    for(d in unique(summary$Var1)){
      num_cell = summary[summary$Var1 == d,]$Freq[1]
      stat_donors[rindex,d] = num_cell
    }
    stat_donors[rindex,1] = row_name
    stat_donors[is.na(stat_donors)] = 0
    assignment[,row_name] = merge(assignment, donorsid, by.x = "barcode", by.y = "cell", all = TRUE)$donor_id
    
    
    f = strsplit(f, split='/vireo+',fixed=TRUE)
    if(length(f[[1]]) == 1){
        f = f[[1]][1]
    }
    else{
	      f = f[[1]][2]
    }
   
    trial_param = strsplit(f, '+',fixed = TRUE)[[1]]
    trial[, paste("Vireo Trial",rindex)] = trial_param 
    
    rindex = rindex + 1
    
  }
  assignment = assignment[order(assignment$barcode),]
  write.table(result_df, paste0(outdir, "vireo_result.csv"), row.names=FALSE, quote = FALSE, sep = ';')
  write.table(stat_donors, paste0(outdir, "vireo_result_detail.csv"), row.names=FALSE, quote = FALSE, sep = ';')
  write.table(assignment, paste0(outdir, "vireo_assignment.csv"), row.names=FALSE, quote = FALSE, sep = ';')
  write.table(trial, paste0(outdir, "vireo_trial.csv"), row.names=FALSE, quote = FALSE, sep = ';')
  write.table(doublet_info, paste0(outdir, "vireo_doublet_info.csv"), row.names=FALSE, quote = FALSE, sep = ';')

}


# analysis souporcell data ----------------------------------------------------------------------------------------------------------------------------------------------------------

num_stat_souporcell <-function(souporcell_list, barcode){
  result_df = as.data.frame(matrix(nrow= length(souporcell_list), ncol=4))
  colnames(result_df) = c("Trial","SNG","DBL","AMB")
  
  f1_cluster = fread(paste(souporcell_list[1],"/clusters.tsv",sep = ""))
  f1_cluster[f1_cluster$status == "doublet",]$assignment = "DBL"
  f1_cluster[f1_cluster$status == "unassigned",]$assignment = "AMB"

  donors = unique(append(f1_cluster$assignment,c("DBL","AMB")))
  
  stat_donors = as.data.frame(matrix(nrow= length(souporcell_list), ncol=length(donors)+1))
  colnames(stat_donors) = append("Trial",donors)
  
  num_barcode = length(barcode)
  assignment = as.data.frame(matrix(ncol= 1, nrow=num_barcode))
  
  assignment[,1] = barcode
  colnames(assignment) = c("barcode")
  
  doublet_info = as.data.frame(matrix(ncol= 1, nrow=num_barcode))
  doublet_info[,1] = barcode
  colnames(doublet_info) = c("barcode")
  
  trial = as.data.frame(matrix(nrow=15, ncol=1))
  colnames(trial) = c("Parameter")
  trial[,1] = c("bam","barcode","fasta","threads","clusters","ploidy","min-alt","min-ref","max-loci","restarts","with-common-variant","with-known-genotype","with-knowngenotypes-sample","skip-remap","ignore")
  
  rindex = 1
  
  for (f in souporcell_list){
    
    cluster = fread(paste(f,"/clusters.tsv",sep = ""))
    
    row_name = paste("Souporcell Trial",rindex)
    
    doublet_info[,c(row_name)] = NA
    doublet_list = cluster[cluster$status%in% c("unassigned","doublet"),c('barcode','assignment')]
    doublet_info[doublet_info$barcode %in% doublet_list$barcode,c(row_name)] = doublet_list[,"assignment"]
    doublet_info[,c(row_name)] = str_replace_all(doublet_info[,c(row_name)],"\\/",",")
    
    cluster[cluster$status == "doublet",]$assignment = "DBL"
    cluster[cluster$status == "unassigned",]$assignment = "AMB"

    num_sin = nrow(cluster[cluster$status=="singlet",])
    num_dou = nrow(cluster[cluster$status=="doublet",])
    num_amb = nrow(cluster[cluster$status=="unassigned",])

    result_df[rindex,1] = row_name
    result_df[rindex,2] = num_sin
    result_df[rindex,3] = num_dou
    result_df[rindex,4] = num_amb
    
    donorlist = unique(append(cluster$assignment,c("DBL","AMB")))
    
    for(d in donorlist){
      num_cell = nrow(cluster[cluster$assignment == d,])
      stat_donors[rindex,d] = num_cell
    }
    stat_donors[rindex,1] = row_name
    stat_donors[is.na(stat_donors)] = 0

    assignment[,row_name] = merge(assignment, cluster, by = "barcode", all = TRUE)$assignment
    
    
    f = strsplit(f, split='/soup+',fixed=TRUE)
    if(length(f[[1]]) == 1){
        f = f[[1]][1]
    }else{
	      f = f[[1]][2]
    }
    
    trial_param = strsplit(f, '+',fixed = TRUE)[[1]]
    trial[, paste("Souporcell Trial",rindex)] = trial_param 
    
    rindex = rindex + 1
    
  }
  assignment = assignment[order(assignment$barcode),]
  stat_donors = stat_donors[, c("Trial", sort(setdiff(names(stat_donors), "Trial")))]

  write.table(result_df, paste0(outdir, "souporcell_result.csv"), row.names=FALSE, quote = FALSE, sep = ';')
  write.table(stat_donors, paste0(outdir, "souporcell_result_detail.csv"), row.names=FALSE, quote = FALSE, sep = ';')
  write.table(assignment, paste0(outdir, "souporcell_assignment.csv"), row.names=FALSE, quote = FALSE, sep = ';')
  write.table(trial, paste0(outdir, "souporcell_trial.csv"), row.names=FALSE, quote = FALSE, sep = ';')
  write.table(doublet_info, paste0(outdir, "souporcell_doublet_info.csv"), row.names=FALSE, quote = FALSE, sep = ';')
  
}


#analysis scSplit data ----------------------------------------------------------------------------------------------------------------------------------------------------------

num_stat_scSplit <-function(scSplit_list, barcode){

  result_df = as.data.frame(matrix(nrow= length(scSplit_list), ncol=4))
  colnames(result_df) = c("Trial","SNG","DBL","AMB")
    
  f1_result = fread(paste(scSplit_list[1],"/scSplit_result.csv",sep = ""))
  f1_prob = fread(paste(scSplit_list[1],"/scSplit_P_s_c.csv",sep = ""))
  
  f1_result = transform(f1_result, new=do.call(rbind, strsplit(Cluster, '-', fixed=TRUE)), stringsAsFactors=F)
  f1_result[,2] = f1_result[,4]
  colnames(f1_result) = c("Barcode",'Donor',"Droplet","Cluster")
  f1_result[f1_result$Droplet == "DBL"]$Donor = "DBL"
  
  donors = unique(append(f1_result$Donor,c("DBL","AMB")))
  stat_donors = as.data.frame(matrix(nrow= length(scSplit_list), ncol=length(donors)+1))
  colnames(stat_donors) = append("Trial",donors)
  
  num_barcode = length(barcode)
  assignment = as.data.frame(matrix(ncol= 1, nrow=num_barcode))
  
  assignment[,1] = barcode
  colnames(assignment) = c("barcode")
  
  trial = as.data.frame(matrix(nrow=8, ncol=1))
  colnames(trial) = c("Paramter")
  trial[,1] = c("vcf","bam", "common-snp","expected-num-samples","max-num-subpopulations-autodetect","num-EM","correction-doublets","known-genotypes")
  
  
  rindex = 1
  for (f in scSplit_list){
    result = fread(paste(f,"/scSplit_result.csv",sep = ""))
    prob = fread(paste(f,"/scSplit_P_s_c.csv",sep = ""))
    
    result = transform(result, new=do.call(rbind, strsplit(Cluster, '-', fixed=TRUE)), stringsAsFactors=F)
    result[,2] = result[,4]
    colnames(result) = c("Barcode",'Donor',"Droplet","Cluster")
    result[result$Droplet == "DBL"]$Donor = "DBL"
  
    num_dou = nrow(result[result$Droplet == "DBL",])
    num_amb = 0
    num_sin = nrow(result[result$Droplet == "SNG",])

    row_name = paste("scSplit Trial",rindex)
    result_df[rindex,1] = row_name
    result_df[rindex,2] = num_sin
    result_df[rindex,3] = num_dou
    result_df[rindex,4] = num_amb
    
    for(d in unique(result$Donor)){
      num_cell = nrow(result[result$Donor == d,])
      stat_donors[rindex,d] = num_cell
    }
    stat_donors[rindex,1] = row_name
    stat_donors[is.na(stat_donors)] = 0
    
    assignment[,row_name] = merge(assignment, result, by.x = "barcode", by.y = "Barcode", all = TRUE)$Donor
    f = strsplit(f, split='/scSplit+',fixed=TRUE)
    if(length(f[[1]]) == 1){
      f = f[[1]][1]
    }else{
      f = f[[1]][2]
    }
    trial_param = strsplit(f, '+',fixed = TRUE)[[1]]
    trial[, paste("scSplit Trial",rindex)] = trial_param 
    
    rindex = rindex + 1
    
  }
  assignment = assignment[order(assignment$barcode),]
  write.table(result_df, paste0(outdir, "scSplit_result_ana.csv"), row.names=FALSE, quote = FALSE, sep = ';')
  write.table(stat_donors, paste0(outdir, "scSplit_result_detail.csv"), row.names=FALSE, quote = FALSE, sep = ';')
  write.table(assignment, paste0(outdir, "scSplit_assignment.csv"), row.names=FALSE, quote = FALSE, sep = ';')
  write.table(trial, paste0(outdir, "scSplit_trial.csv"), row.names=FALSE, quote = FALSE, sep = ';')
  
}
 
parser <- ArgumentParser()
parser$add_argument("--file", nargs=1, help="demultiplexing output file")
parser$add_argument("--barcode", nargs=1, help="barcode file")
parser$add_argument("--outdir", nargs=1, help="Output Directory")
args <- parser$parse_args()
file <- args$file
barcode <- fread(args$barcode, header = FALSE)
barcode <- barcode$V1
outdir <- args$outdir

file <- strsplit(file, ':',fixed = TRUE)[[1]]
num_tool = 0
demux_list= str_subset(file,'demux')
if (length(demux_list) != 0){
  num_stat_demux(demux_list, barcode)
  num_tool = num_tool + 1
}else{
  print("No demuxlet output found!")
}

freemuxlet_list= str_subset(file,'freemuxlet')

if (length(freemuxlet_list) != 0){
  freemuxlet_list= paste0(freemuxlet_list, "/freemuxlet.clust1.samples.gz")
  num_stat_freemuxlet(freemuxlet_list, barcode)
  num_tool = num_tool + 1
}else{
  print("No freemuxlet output found!")
}

vireo_list = str_subset(file,'vireo')
if (length(vireo_list) != 0){
  num_stat_vireo(vireo_list, barcode)
  num_tool = num_tool + 1
}else{
  print("No vireo output found!")
}
  
souporcell_list = str_subset(file,'soup')
if (length(souporcell_list) != 0){
  num_stat_souporcell(souporcell_list, barcode)
  num_tool = num_tool + 1
}else{
  print("No souporcell output found!")
}

scSplit_list = str_subset(file,'scSplit')
if (length(scSplit_list) != 0){
  num_stat_scSplit(scSplit_list, barcode)
  num_tool = num_tool + 1
}else{
  print("No scSplit output found!")
}

#assignment_all <- list.files(path =  ".", pattern = "_assignment.csv", full.names = TRUE) %>%
assignment_all <- list.files(path = outdir, pattern = "_assignment.csv", full.names = TRUE) %>%
    lapply(fread) %>%
    bind_cols()  %>% 
    as.data.frame()
colnames(assignment_all)[1] = "Barcode"
if (num_tool > 1){
  assignment_all = assignment_all[,-which(startsWith(colnames(assignment_all),'barcode'))]
}
write.csv(assignment_all,paste0(outdir, "assignments.csv"))
