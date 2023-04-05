# -------------------------- #
# Prepare data for stacked barplot
# -------------------------- #

#' @name
#' GenomewideStackedBar
#'
#' @title
#' Prep data for genome wide stacked barplot
#' Create stacked bar plot
#'
#' @description
#' Take output from bedtools intersect and prepare data for stacked barplot
#' Create stacked barplot for 6 different CNA types
#' homdel, loh, otherloss, nochange, gain, biggain
#'
#' @param
#' filestub = dir (with trailling slash) containing *_partitioned_regions_all_CNA_types_overlapped.out files for all cna types
#' segfile_name = label of cohort, eg 'TGCT'
#'
#' @return
#' Prepped data
#' Stacked barplot
#'

# collate all segments in subclones files across cohort
# BA edited function below
GenomewideStackedBar <- function (segfile_dir, segfile_name) 
{
  # load libraries
  library(reshape2)

  all <- read.table(file.path(segfile_dir, paste0(segfile_name, "_partitioned_regions_all_CNA_types_overlapped.out")), 
    sep = "\t", stringsAsFactors = F)
  
  data <- all[, c(1,2,3,8,9)]
  
  colnames(data) <- c("poschr", "posleft", "posright", "value", "region")
  
  data$value <- as.numeric(data$value)
  
  data_wide <- dcast(data, poschr + posleft + posright ~ region, value.var = "value")
  
  data_wide$chr_num <- as.numeric(gsub("chr", "", data_wide$poschr))
  
  data_wide$posleft <- as.numeric(data_wide$posleft)
  
  data_wide$posright <- as.numeric(data_wide$posright)
  
  data_wide <- data_wide[order(data_wide$chr_num, data_wide$posleft, data_wide$posright), ]
  
  # do we need row.names?
  row.names(data_wide) <- paste(data_wide$poschr, data_wide$posleft, data_wide$posright, sep = "_")
  
  data_wide <- data_wide[, c("poschr", "posleft", "posright", "homdel", "loh", "otherloss", "nochange", "gain", "biggain")]
  
  pos <- paste(all[, 1], all[, 2], all[, 3], sep = "_")
  all <- cbind(all, pos)
  pos <- unique(pos)
  
  poschr <- sapply(pos, function(x) {strsplit(x, "_")[[1]][1]})
  
  posleft <- sapply(pos, function(x) {strsplit(x, "_")[[1]][2]})
  
  posright <- sapply(pos, function(x) {strsplit(x, "_")[[1]][3]})
  
  forplot <- cbind(poschr, posleft, posright)
  homdel <- rep(0, nrow(forplot))
  loh <- rep(0, nrow(forplot))
  otherloss <- rep(0, nrow(forplot))
  nochange <- rep(0, nrow(forplot))
  gain <- rep(0, nrow(forplot))
  biggain <- rep(0, nrow(forplot))
  
  forplot <- as.data.frame(cbind(forplot, homdel, loh, otherloss, 
                nochange, gain, biggain))
  
  for (i in 2:ncol(forplot)) {
    forplot[, i] <- as.integer(as.character(forplot[, i]))
  }
  
  cnas <- c("homdel", "loh", "otherloss", "nochange", "gain", "biggain")
  
  for (cna in 1:length(cnas)) {
    allcna <- all[which(all[, 9] == cnas[cna] & all[, 8] != 0), ]
    test <- match(rownames(forplot), allcna$pos)
    forplot_rows <- which(!is.na(test))
    allcna_rows <- test[-which(is.na(test))]
    
    forplot[forplot_rows, cnas[cna]] <- allcna[allcna_rows, 8]
  }
  
  write.table(data_wide, file.path(segfile_dir, paste0(segfile_name, "_partitioned_regions_all_CNA_types_overlapped_forplot_alt.out")), 
          quote = F, sep = "\t")
  
  write.table(forplot, file.path(segfile_dir, paste0(segfile_name, "_partitioned_regions_all_CNA_types_overlapped_forplot.out")), 
          quote = F, sep = "\t")
}
