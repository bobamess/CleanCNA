# -------------------------- #
# Prepare data for stacked barplot
# -------------------------- #

#' @name
#' GenomewideStackedBarPlot
#'
#' @title
#' Create stacked bar plot
#'
#' @description
#' Create stacked barplot for 6 different CNA types
#' homdel, loh, otherloss, nochange, gain, biggain, split into loss, no change, and gain
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
GenomewideStackedBarPlot <- function (segfile_dir, segfile_name, number_samples) 
{
  # load libraries
  library(scales)
  library(reshape2)
  library(ggplot2)
  library(tools) # for toTitleCase()
  library(grDevices) # for png()
  # png() also requires options(device='x11', bitmapType="cairo") to access x11
  # - this can be put in .Rprofile and enabled with the following in .bash_profile
  # R_PROFILE_USER=/users/$USER/.Rprofile
  # export R_LIBS_USER
  
  number_samples <- as.numeric(number_samples)
  
  toplot <- read.table(file.path(segfile_dir, paste0(segfile_name, "_partitioned_regions_all_CNA_types_overlapped_forplot.out")), 
      header = T, stringsAsFactors = F, fill = T)
      
  centro <- data.frame(pos = c(123400000, 93900000, 90900000, 
          5e+07, 48750000, 60550000, 60100000, 45200000, 43850000, 
          39800000, 53400000, 35500000, 17700000, 17150000, 1.9e+07, 
          36850000, 25050000, 18450000, 26150000, 28050000, 11950000, 
          15550000), poschr = c(1:22))
              
  colours <- c("#FF0000FF", "#FF9E81FF", "wheat1", "#A3A9FFFF", "#7863FFFF", "#0000FFFF") # to match heatmap including "nochange" which is no longer white
  
  biggaincol <- colours[1]
  gaincol <- colours[2]
  nochangecol <- colours[3]
  otherlosscol <- colours[4]
  lohcol <- colours[5]
  homdelcol <- colours[6]
  toplot$valuenochange <- toplot$nochange
  toplot$valuegain <- toplot$gain
  toplot$valuegain2 <- toplot$gain + toplot$biggain

  # BA changed order to match order in legend colour scale
  toplot$valueloss <- toplot$homdel
  toplot$valueloss2 <- toplot$homdel + toplot$loh
  toplot$valueloss3 <- toplot$homdel + toplot$loh + toplot$otherloss
  
  # convert counts to percentages
  
  for (column in 4:15) {
    toplot[, column] <- (toplot[, column]/number_samples) * 100
  }
  
  toplot$poschr <- as.integer(gsub("chr", "", toplot$poschr))
  
  # next for loop seems redundant?
  all <- toplot[which(toplot$poschr == 1), ]
  
  for (chr in 2:23) {
    sub <- toplot[which(toplot$poschr == chr), ]
    all <- rbind(all, sub)
  }
  
  toplot <- all

  gainstoplot <- toplot[, c("poschr", "posleft", "posright", "gain", "biggain", "valuegain", "valuegain2")]
      
  gains <- ggplot(gainstoplot) + 
    geom_rect(aes(xmin = posleft, xmax = posright, ymin = 0, ymax = valuegain, fill = gaincol)) + 
    geom_rect(aes(xmin = posleft, xmax = posright, ymin = valuegain, ymax = valuegain2, fill = biggaincol)) + 
    ylab("% tumour samples") + 
    facet_grid(~poschr, scales = "free_x", space = "free_x") + 
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
    coord_cartesian(ylim = c(0, 100)) + 
    scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100), labels = c(0, 20, 40, 60, 80, 100)) + 
    scale_x_continuous(breaks = NULL) + 
    geom_point(data = centro, aes(x = pos, y = 0), size = 1.5) + 
    geom_hline(yintercept = 0, colour = "black", size = 0.5) + 
    scale_fill_manual(values = c(biggaincol, gaincol), labels = c("big gain", "gain")) + 
    guides(fill = guide_legend(title = "Gain"))
    
  pdf(file.path(segfile_dir, paste0(segfile_name, "_genomewide_stackedbar_gains.pdf")), width = 10, height = 2)
    print(gains)
  dev.off()
  
  options(device='x11', bitmapType="cairo")
  
  png(file.path(segfile_dir, paste0(toTitleCase(segfile_name), "_gains_mqc.png")), width = 960, height = 192)    
    print(gains)
  dev.off()
  
  lossestoplot <- toplot[, c("poschr", "posleft", "posright", "homdel", "loh", "otherloss", "valueloss", "valueloss2", "valueloss3")]
  
  alt_lossestoplot <- melt(toplot[, c("poschr", "posleft", "posright", "homdel", "loh", "otherloss")], variable.name = "Losses", id.var = c("poschr", "posleft", "posright"), value.name = "pct")
  
  alt_lossestoplot$pct <- as.numeric(alt_lossestoplot$pct)
  
  alt_lossestoplot$Losses <- as.factor(alt_lossestoplot$Losses)
  
  # see https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin  
  # see also https://ggplot2.tidyverse.org/reference/scale_manual.html
  
  losses_col <- c("otherlosscol" = otherlosscol, "lohcol" = lohcol, "homdelcol" = homdelcol)
  
  Losses <- alt_lossestoplot$Losses
  
  losses <- ggplot(lossestoplot) + 
          # BA changed below to match order in legend colour scale
          geom_rect(aes(xmin = posleft, xmax = posright, ymin = 0, ymax = valueloss, fill = homdelcol)) + 
          geom_rect(aes(xmin = posleft, xmax = posright, ymin = valueloss, ymax = valueloss2, fill = lohcol)) + 
          geom_rect(aes(xmin = posleft, xmax = posright, ymin = valueloss2, ymax = valueloss3, fill = otherlosscol)) + 
          
          ylab("% tumour samples") + 
          facet_grid(~poschr, scales = "free_x", space = "free_x") + 
          theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
          coord_cartesian(ylim = c(0, 100)) + 
          scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100), labels = c(0, 20, 40, 60, 80, 100)) + 
          scale_x_continuous(breaks = NULL) + 
          geom_point(data = centro, aes(x = pos, y = 0), size = 1.5) + 
          geom_hline(yintercept = 0, colour = "black", size = 0.5) + 
          scale_fill_manual(values = c(homdelcol, lohcol, otherlosscol), labels = c("homdel", "loh", "other loss")) +       
          guides(fill = guide_legend(title = "Loss"))
          
  # reversing the order of values and labels in scale_fill_manual, not only reverses the legend order but also corrects the colours in the plot to match those in heatmap, 
  # but why should this be as colours for the plot should be defined by fill attribute, which were already correct? - seems to be a bug - is it to do with levels?
    
  alt_losses <- ggplot(alt_lossestoplot, aes(x = posleft, y = pct)) +  
          
          geom_bar(position="stack", stat = "identity") + 
          
          scale_colour_manual(values = losses_col) + 
          
          ylab("% tumour samples") + 
          facet_grid(~poschr, scales = "free_x", space = "free_x") + 
          theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
          coord_cartesian(ylim = c(0, 100)) + 
          scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100), labels = c(0, 20, 40, 60, 80, 100)) + 
          scale_x_continuous(breaks = NULL) + 
          geom_point(data = centro, aes(x = pos, y = 0), size = 1.5) + 
          geom_hline(yintercept = 0, colour = "black", size = 0.5) +     
          guides(fill = guide_legend(title = "Loss"))
          

  pdf(file.path(segfile_dir, paste0(segfile_name, "_genomewide_stackedbar_losses.pdf")), width = 10, height = 2)
    print(losses)
  dev.off()
  
  png(file.path(segfile_dir, paste0(toTitleCase(segfile_name), "_losses_mqc.png")), width = 960, height = 192)
    print(losses)
  dev.off()
  
  nochangetoplot <- toplot[, c("poschr", "posleft", "posright", "nochange", "valuenochange")]
  
  nochangetoplot$value = nochangetoplot$nochange # nochangetoplot$value not used?
  
  nochange <- ggplot(nochangetoplot) + 
        geom_rect(aes(ymin = 0, xmin = posleft, xmax = posright, ymax = valuenochange, fill = nochangecol)) + 
        ylab("% tumour samples") + 
        facet_grid(~poschr, scales = "free_x", space = "free_x") + 
        theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
        coord_cartesian(ylim = c(0, 100)) + 
        scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100), labels = c(0, 20, 40, 60, 80, 100)) + 
        scale_x_continuous(breaks = NULL) + 
        geom_point(data = centro, aes(x = pos, y = 0), size = 1.5) + 
        geom_hline(yintercept = 0, colour = "black", size = 0.5) + 
        scale_fill_manual(values = nochangecol, labels = "no change") + 
        guides(fill = guide_legend(title = NULL))
  
  pdf(file.path(segfile_dir, paste0(segfile_name, "_genomewide_stackedbar_nochange.pdf")), width = 10, height = 2)
    print(nochange)
  dev.off()
  
  png(file.path(segfile_dir, paste0(toTitleCase(segfile_name), "_no_change_mqc.png")), width = 960, height = 192)
    print(nochange)
  dev.off()
}

