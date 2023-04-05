# -------------------------- #
# Label bedtools_genomecov files with type of CNA
# -------------------------- #

#' @name
#' LabelGenomecov
#'
#' @title
#' Collate all segments in all samples in the cohort
#'
#' @description
#' Label bedtools_genomecov files with type of CNA
#'
#' @param
#' filestub = dir (with trailling slash) containing *regions_normalised_out files for all cna types
#' segfile_name = label of cohort, eg 'TGCT'
#'
#' @return
#' writes same name file with cna types labelled
#'

# collate all segments in subclones files across cohort
LabelGenomecov <- function (filestub, segfile_name) 
{
	cnas = c("homdel", "loh", "otherloss", "nochange", "gain", "biggain")
	
	for (i in 1:length(cnas)) {
		# BA edited line below
		file = read.table(file.path(filestub, paste0(segfile_name, "_", cnas[i], "_regions.out")), sep = "\t", stringsAsFactors = F)
		file = cbind(file, cnas[i])
		
		# BA edited line below
		write.table(file, file.path(filestub, paste0(segfile_name, "_", cnas[i], "_regions.out")), sep = "\t", col.names = F, row.names = F, quote = F)
	}
}
