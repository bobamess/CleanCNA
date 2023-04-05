# -------------------------- #
# Prepare data for stacked barplot
# -------------------------- #

#' @name
#' CleanBedopsOutput
#'
#' @title
#' Clean bedops output
#'
#' @description
#' Cleans output from bedops which leaves segments spanning whole chrs which screw up downstream steps
#'
#' @param
#' filestub = dir (with trailling slash) containing *_bedops
#' segfile_name = label of cohort, eg 'TGCT'
#'
#' @return
#' Clean version of bedops file
#'

# collate all segments in subclones files across cohort
# BA edited function below replacing filestub with segfile_dir
CleanBedopsOutput <- function (segfile_dir, segfile_name)
{
  bedops = read.csv(file.path(segfile_dir, paste0(segfile_name, "_bedops")), 
            sep = "\t", header = F, stringsAsFactors = F)
            
  bedops[, 1] = gsub("chr", "", bedops[, 1])
  bedops[, 4] = paste(bedops[, 1], bedops[, 2], bedops[, 3], sep = "_")
  
  for (i in 1:3) {
    bedops[, i] = as.integer(bedops[, i])
  }
  
  toremove = which(duplicated(bedops))
  
  if (length(toremove) > 0) {
    bedops = bedops[-toremove, ]
  }
  
  lastrows = c()
  
  for (chr in 1:23) {
    thischr = bedops[which(bedops[, 1] == chr), ]
    rightmost = which(thischr[, 2] == 0 & thischr[, 3] >= max(thischr[, 3]))
    
    if (length(rightmost) > 0) {
      print(chr)
      lastrows[chr] = which(bedops[, 4] == thischr[rightmost, 4])
    } else if (length(rightmost) == 0) {
      print(chr)
      lastrows[chr] = 0
    }
  }
  
  lastrows = lastrows[-which(lastrows == 0)]
  
  if (length(lastrows) > 0) {
    bedops = bedops[-lastrows, ]
  }
  
  bedops = bedops[order(bedops[, 1], bedops[, 2], bedops[, 3]), ]
  bedops[, 1] = paste("chr", bedops[, 1], sep = "")
  
  write.table(bedops[, 1:3], file.path(segfile_dir, paste0(segfile_name, "_bedops")), 
        sep = "\t", row.names = F, quote = F, col.names = F)
}
