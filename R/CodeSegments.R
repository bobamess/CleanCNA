# -------------------------- #
# Code segments from subclones files
# -------------------------- #

#' @name
#' CodeSegments
#'
#' @title
#' Code all segments in all samples in the cohort
#'
#' @description
#' Codes all segments as a particular CNA type from across the cohort
#'
#' @param
#' segfile_dir = where segfile will go
#' segfile_name = cohort name
#'
#' @return
#' A table for use with bedtools etc
#'
#'



# code all segments into homdel, loh, otherloss, no change, gain, big gain
CodeSegments <- function (segfile_dir, segfile_name) 
{
  # BA changed line below
  subs <- read.table(file.path(segfile_dir, paste0(segfile_name, "_segsfull.txt")), sep = "\t", header = T, stringsAsFactors = F)
  
  subs$total_cn = rowSums(subs[c("nMajor", "nMinor")])
  subs$class[subs$class == "dip"] <- 2
  subs$class[subs$class == "tetra"] <- 4
  colnames(subs)[14] = "dip.tetra"
  subs$dip.tetra = as.numeric(subs$dip.tetra)
  subs$coded_total_cn = NA
  subs$coded_total_cn[which(subs$total_cn == 0)] = "homdel"
  subs$coded_total_cn[which(subs$nMajor > 0 & subs$nMinor == 0)] = "loh"
  
  subs$coded_total_cn[which(subs$total_cn < subs$dip.tetra & 
    subs$total_cn != 0 & subs$nMinor != 0)] = "otherloss"
  
  subs$coded_total_cn[which((subs$nMajor == subs$nMinor | subs$nMajor == 3 & subs$nMinor == 1) & (subs$total_cn == subs$dip.tetra))] = "nochange"
  
  subs$coded_total_cn[which(subs$total_cn > subs$dip.tetra)] = "gain"
  subs$coded_total_cn[which(subs$total_cn > (5 * subs$dip.tetra))] = "biggain"
  subs$chr[subs$chr == "X"] = 23
  subs <- subs[order(subs$chr, subs$startpos, subs$endpos), ]
  subs$chr = paste0("chr", subs$chr)
  
  cnv_levels <- c("homdel", "loh", "otherloss", "nochange", "gain", "biggain")
  
  cnv_numbers = c(-2, -1, -0.5, 2, 5, 10)
  for (cnv in 1:length(cnv_levels)) {
    cnv_subs = subs[which(subs$coded_total_cn == cnv_levels[cnv]), ]
    cnv_subs = cnv_subs[order(cnv_subs$chr, cnv_subs$startpos, cnv_subs$endpos), ]
    
    # BA edited line below
    write.table(cnv_subs, file.path(segfile_dir, paste0(segfile_name, "_", cnv_levels[cnv], ".full")), quote = F, row.names = F, sep = "\t")
  
    subs_bed = cbind(cnv_subs[c("chr", "startpos", "endpos", "sample")])
    colnames(subs_bed)[4] = "name"
    print(nrow(subs_bed))
    
    # BA edited line below
    write.table(subs_bed, file.path(segfile_dir, paste0(segfile_name, "_", cnv_levels[cnv], ".bed")), quote = F, row.names = F, col.names = F, sep = "\t")
  }
}
