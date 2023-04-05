# -------------------------- #
# Make heatmap
# -------------------------- #

#' @name
#' CNAsDriversHeatmapDendrogram
#'
#' @title
#' Heatmap to show CNA status in a set of driver genes across a cohort of tumour samples
#' Dendrogram to go alongside - heatmap ordered to match dendrogram
#'
#' @description
#' Create heatmap showing CNA status at driver genes using ggplot2
#' homdel, loh, otherloss, nochange, gain, biggain
#'
#' @param
#' segfile_dir = dir with all *_subclones.txt files are
#' segfile_name = label of cohort, eg 'TGCT'
#' bed_file = tab sep file, 7 cols, chr,start,end,strand,ens,ens,genename
#' driver_file = tab sep file, colnames SYMBOL (of gene) and TIER (1,2,etc) existing
#'
#' @return
#' Heatmap and associated dendrogram (in same order) of CNAs v drivers
#'

# collate all segments in subclones files across cohort
# BA edited function below
CNAsDriversHeatmapDendrogram <- function (segfile_dir, segfile_name, bed_file, driver_file, metadata_file)
{
  # load libraries
  library(dplyr)
  library(ggplot2)
  library(ggdendro)
  library(reshape2)
  library(grid)
   
  subclones_dir <- segfile_dir
  
  subclones <- system(paste0("ls ", file.path(subclones_dir, "*subclones.txt")), intern = T)
  
  samples <- basename(subclones)
  
  samples <- gsub("_subclones.txt", "", samples)
  
  bed <- read.table(bed_file, stringsAsFactors = F)
  
  dri <- read.table(driver_file, header = T, stringsAsFactors = F, fill = T)
   
  # BA changed line below to accommodate alternative TIER values
  
  dri <- dri[which(as.character(dri$TIER) %in% c("1", "2", "3", "High", "Mid", "Low")), ]
  
  bed <- merge(bed, dri, by.x = c(7), by.y = c("SYMBOL")) # this changes order of columns
  
  bed[which(bed[,2] == "chrX"), 2] <- "chr23"
  
  bed[,1] <- paste(bed[,2], bed[,1], sep = " ")
  
  # remove chromosomes/scaffolds not in data
  bed[, 9] <- as.integer(gsub("chr", "", bed[, 2]))
  
  bed <- bed[which(!is.na(bed[, 9])), c(1:8)] 

  # concatenate chromosome, gene_name (these 2 are already concat in col 1) and transcript_id
  bed[,9] <- paste(bed[,1], bed[,7], sep = " ") 
  
  positions_chr <- as.character(bed[, 2])
  positions_chr <- as.character(gsub("chr", "", positions_chr))
  positions_start <- as.character(bed[, 3])
  positions_end <- as.character(bed[, 4])
  positions_gene <- as.character(bed[, 9]) # actually this is now a transcript_id in bed_file
  
  genes <- as.data.frame(cbind(positions_gene, positions_chr, positions_start, positions_end))
  
  for (i in 1:ncol(genes)) {
  genes[, i] = as.character(genes[, i])
  }
  
  genes <- genes[order(as.integer(genes$positions_chr)), ]
  genes$positions_chr <- as.integer(genes$positions_chr)
  genes$positions_start <- as.integer(genes$positions_start)
  genes$positions_end <- as.integer(genes$positions_end)
  
  genes <- genes[which(!is.na(genes[,c("positions_chr")])),]
  
  subs <- read.table(file.path(segfile_dir, paste0(segfile_name, "_segsfull.txt")), sep = "\t", stringsAsFactors = F, hea = T)
  
  subs$class[subs$class == "dip"] <- 2
  subs$class[subs$class == "tetra"] <- 4
  subs$class <- as.integer(subs$class)
  subs$coded_total_cn <- NA
  subs$total_cn <- subs$nMajor + subs$nMinor
  subs$coded_total_cn[which(subs$total_cn == 0)] <- "homdel"
  subs$coded_total_cn[which(subs$nMajor > 0 & subs$nMinor == 0)] <- "loh"
  subs$coded_total_cn[which(subs$total_cn < subs$class & subs$total_cn != 0 & subs$nMinor != 0)] <- "otherloss"
  subs$coded_total_cn[which((subs$nMajor == subs$nMinor | subs$nMajor == 3 & subs$nMinor == 1) & (subs$total_cn == subs$class))] <- "nochange"
  subs$coded_total_cn[which(subs$total_cn > subs$class)] <- "gain"
  subs$coded_total_cn[which(subs$total_cn > (5 * subs$class))] <- "biggain"
  subs$chr[subs$chr == "X"] <- 23
  
  # in lines below original code used name tomelt instead of mat
  
  mat <- matrix(NA, nrow = length(unique(subs$sample)), ncol = nrow(genes)) # create matrix with NA values with rows for each sample and cols for each gene (actually transcript_id)
  rownames(mat) <- unique(subs$sample)
  colnames(mat) <- genes$positions_gene
  
  subs$chr <- as.integer(subs$chr)

  for (i in 1:length(samples)) {
    sub <- subs[subs$sample == samples[i], ] # subclones for sample[i]
    
    for (g in 1:nrow(genes)) {
      chr <- sub[sub$chr == genes$positions_chr[g], ] # subclones for sample[i] with same chromosome as genes[g] (actually transcript_id), if any
      
      if (nrow(chr) == 0) { # if none, go to next gene
        next
      }
      
      row <- which((chr$startpos <= genes$positions_start[g]) & # find row(s) for subclones selected that have both startpos and endpos equal to or beyond gene (actually transcript_id)
        (chr$endpos >= genes$positions_end[g])) 
      
      if (length(row) == 0) {
        minrows <- max(which(chr$startpos <= genes$positions_start[g]), na.rm = TRUE) # find row in chr that has startpos closest to gene positions_start, if any 
        maxrows <- min(which(chr$endpos >= genes$positions_end[g]), na.rm = TRUE) # find row in chr that has endpos closest to gene positions_end, if any 
        if (minrows == -Inf & nrow(chr) == 1) {
          minrows <- 1
        }
        # how do you get -Inf or Inf?
        if (maxrows == Inf & nrow(chr) == 1) {
          maxrows <- 1
        }
        rows <- minrows:maxrows
        longestregion <- which(chr[rows, c("endpos")] - chr[rows, c("startpos")] == max(chr[rows, c("endpos")] - chr[rows, c("startpos")], na.rm = TRUE))
        row <- rows[longestregion]
        
        mat[i, g] <- chr$coded_total_cn[row]
      } else {
        mat[i, g] <- chr$coded_total_cn[row]
      }
    }
  }
  
  mat[mat == "homdel"] <- 0
  mat[mat == "loh"] <- 1
  mat[mat == "otherloss"] <- 2
  mat[mat == "nochange"] <- 3
  mat[mat == "gain"] <- 4
  mat[mat == "biggain"] <- 5
  
  for (i in 1:ncol(mat)) {
    mat[, i] <- as.numeric(mat[, i])
  }
  
  write.table(mat, file.path(segfile_dir, paste0(segfile_name, "_dendrogram.csv")), 
    sep = "\t", quote = F)
  
  melted_mat <- melt(mat)
  colnames(melted_mat) <- c("tumour", "driver", "normalised_total_cn")
  
  write.table(melted_mat, file.path(segfile_dir, paste0(segfile_name, "_heatmap.csv")), 
    row.names = F, quote = F)
  
  write.table(mat, file.path(segfile_dir, paste0(segfile_name, "_dendrogram.csv")), 
    sep = "\t", quote = F)
}
