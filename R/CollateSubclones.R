# -------------------------- #
# Collate all segments from subclones files
# -------------------------- #

#' @name
#' CollateSubclones
#'
#' @title
#' Collate all segments in all samples in the cohort
#'
#' @description
#' Prepares subclones info to label copy number of segments
#'
#' @param
#' qc table from PrepDataForSummary
#' segfile_dir = where segfile will go
#' segfile_name = cohort name
#'
#' @return
#' A table with all segments from all subclones files collated
#'
#'

# collate all segments in subclones files across cohort
CollateSubclones <- function (qc, segfile_dir, segfile_name) 
{
  subsfull = matrix(nrow = 0, ncol = 13)
  
  qc <- read.csv(qc, sep = "\t", stringsAsFactors = F)
  
  # BA - below is another way to make a key using "." instead of ":" - different to file names and paths
  samples <- paste0(qc$participant_id, ".", qc$tumour_sample_platekey, ".", qc$germline_sample_platekey)
  
  # BA changed line below
  subclonesPaths <- list.files(segfile_dir, , pattern = "subclones.txt", full.names = T)
  
  # BA changed line below
  purityPaths <- gsub(pattern = "subclones.txt", replacement = "cellularity_ploidy.txt", subclonesPaths)
  
  for (i in 1:length(subclonesPaths)) {
    # BA edited line below
    sub <- read.csv(subclonesPaths[i], header = T, stringsAsFactors = F, sep = "\t")[, c(1:13)]
    
    samplename <- basename(subclonesPaths[i])
    
    samplename <- gsub(pattern = "_subclones.txt", replacement = "", samplename)
    
    sub$subclonal <- sub$frac1_A > sub$frac2_A
    firstsubclone <- sub[, 8:10]
    secondsubclone <- sub[, 11:13]
    sub[which(sub$subclonal == F), 8:10] <- secondsubclone[which(sub$subclonal == F), ]
    sub[which(sub$subclonal == F), 11:13] <- firstsubclone[which(sub$subclonal == F), ]
    sub <- sub[, -c(4:7, 14)]
    
    #sub <- cbind(as.character(samples[i]), sub) # assumes same order for samples as in subclonesPaths, which may not be true!!!
    # instead use samplename as derived from subclonesPaths[i]
    sub <- cbind(as.character(samplename), sub)
    
    colnames(sub)[1:6] <- c("sample", "chr", "startpos", "endpos", "nMajor", "nMinor")
    
    # BA edited line below
    sub$ploidy <- signif(read.table(purityPaths[i], header = T, stringsAsFactors = F, sep = "\t")[1, 3])
    
    sub$size <- sub$endpos - sub$startpos + 1
    sub2 <- sub
    sub2[is.na(sub2)] <- 0
    genome.length <- 2923364639
    sub2fracLOH <- sub2[which(round((sub2$nMinor * sub2$frac1_A) + (sub2$nMin2_A * sub2$frac2_A)) == 0), ]
    sub2fracLOH <- sum(sub2fracLOH$endpos - sub2fracLOH$startpos)/genome.length
    sub$psi_t <- sum(sub2$size * (sub2$nMajor + sub2$nMinor))/sum(sub2$size)
    psi_t <- sub$psi_t
    sampleploidy <- ifelse(psi_t > (sub2fracLOH * -2) + 2.9, "tetra", "dip")
    sub$class <- sampleploidy
    subsfull <- rbind(subsfull, sub)
  }
  
  # BA edited line below
  write.table(subsfull, file.path(segfile_dir, paste0(segfile_name, "_segsfull.txt")), sep = "\t", row.names = F, quote = F)
}
