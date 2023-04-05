# -------------------------- #
# Make complexheatmap
# -------------------------- #

#' @name
#' CNAsDriversComplexHeatmap
#'
#' @title
#' Heatmap to show CNA status in a set of driver genes across a cohort of tumour samples
#' Dendrogram automatically goes alongside
#' Bars classifying samples in lots of ways produced alongside heatmap
#'
#' @description
#' Create heatmap showing CNA status at driver genes using complexheatmap
#' homdel, loh, otherloss, nochange, gain, biggain
#' with various bars down the side
#' for bars down side, feed vectors to function
#'
#' @param
#' segfile_dir = dir with all *_subclones.txt files are
#' segfile_name = label of cohort, eg 'TGCT'
#' plot_name = name for plot
#' tomelt = saved from CNAsDriversHeatmap function
#' bed_file =
#' driver_file = tab sep file, colnames SYMBOL (of gene) and TIER (1,2,etc) existing
#'
#' @return
#' Heatmap, dendrogram of CNAs v drivers, with classifying bars
#'

# BA edited function below
CNAsDriversComplexHeatmap <- function (segfile_dir, segfile_name, plot_name, bed_file, driver_file, metadata_file, sample_list_file, clinical_data_file, last_run_file)
{

  ## load libaries
  library(ComplexHeatmap)
  library(tools) # for toTitleCase()
  library(openxlsx)
  library(circlize)
  library(grDevices) # for png() and colorRampPalette()
  # png() also requires options(device='x11', bitmapType="cairo") to access x11
  # - this can be put in .Rprofile and enabled with the following in .bash_profile
  # R_PROFILE_USER=/users/$USER/.Rprofile
  # export R_LIBS_USER
  
  
  ## read in metadata
  metadata <- read.table(metadata_file, header = T, sep = "\t", stringsAsFactors = F, fill = T)
  
  metadata <- unique(metadata[, c("Source.Label", "Sample", "Time.Point", "Sample.Origin", "Replicate", "Disease", "Group", "Gender", "Age")])
  metadata$Label <- paste0(metadata$Source.Label, "_", metadata$Time.Point, "_", metadata$Sample.Origin, "_",  metadata$Replicate)
  metadata$Label <- gsub("_NA", "", metadata$Label)  # cleans sample names (as metadata$Label) for samples that did not have replicates
  
  row.names(metadata) <- metadata[, "Label"]  
  
  
  ## read in sample list
  sample_list_df <- read.table(sample_list_file, header = T, sep = "\t", stringsAsFactors = F, fill = T)
  
  sample_list_df <- sample_list_df[, c("tumour_sample_platekey", "germline_sample_platekey")]
  
  sample_list_df$tumour_sample_platekey <- gsub(":", "_", sample_list_df$tumour_sample_platekey)
  
  sample_list_df$Reference <- sapply(sample_list_df$germline_sample_platekey, function(x) {unlist(strsplit(x, split = "_"))[3]})
  
  
  ## merge metadata with tumour_sample_platekey 
  metadata <- merge(metadata, sample_list_df, by.x = c("Label"), by.y = c("tumour_sample_platekey"))
  
  
  ## section below repeats what is done in CleanCNA:::CNAsDriversHeatmapDendrogram()
  # so should modify to input directly

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
  
  # in lines below original code used name tomelt instead of mat
  
  mat <- as.matrix(read.table(file.path(segfile_dir, paste0(segfile_name, "_dendrogram.csv")), 
            stringsAsFactors = F, sep = "\t", check.names = FALSE)) # check.names = FALSE retains original column names and 
            # prevents read.table() invoking make.names() which would replace "-" with ".".  This also keeps spaces
  
  rownames(mat) <- gsub(":", "_", rownames(mat)) # may be redundant after changes to init_run.nf
  
  samples <- rownames(mat)
  
  metadata <- metadata[which(metadata$Label %in% samples),]
  
  
  ## read in clinical data
  
  clinical_data_raw <- openxlsx::read.xlsx(clinical_data_file, sheet = "Summary_Metrics")
  
  clinical_data_df <- clinical_data_raw[, c("Patient", "Response_binary", "Histology")]
  
  colnames(clinical_data_df) <- c("Source.Label", "Response", "Histology")
  
  # merge metadata with clinical data
  
  metadata <- merge(metadata, clinical_data_df)
  
  
  ## merge in QC data
  
  qc_data_raw <- read.table(last_run_file, header = T, stringsAsFactors = F, fill = T)
  
  qc_data <- qc_data_raw[, c("tumour_sample_platekey", "QC")]
  
  qc_data$tumour_sample_platekey <- gsub(":", "_", qc_data$tumour_sample_platekey) # may be redundant after changes to init_run.nf
  
  metadata <- merge(metadata, qc_data, by.x = c("Label"), by.y = c("tumour_sample_platekey"))
  
  
  ## how to collapse transcripts for same gene with identical CN profiles across samples to one column in mat?
  
  # transpose mat
  
  mat_t_df <- as.data.frame(t(mat))
  
  # add column with chr gene
  
  mat_t_df$chr_gene <- sapply(rownames(mat_t_df), function(x) {parts <- unlist(strsplit(x, " ", fixed = TRUE)); paste(parts[1], parts[2], sep = ".")})
  
  mat_t_df$gene <- sapply(rownames(mat_t_df), function(x) {parts <- unlist(strsplit(x, " ", fixed = TRUE)); parts[2]})
  
  mat_t_df <- unique(mat_t_df)
  
  mat_t_df$count <- sapply(mat_t_df$chr_gene, function(x) {sum(mat_t_df$chr_gene == x)})
  
  mat_t_df$name <- mat_t_df$chr_gene
  
  mat_t_df$name[which(mat_t_df$count > 1)] <- rownames(mat_t_df)[which(mat_t_df$count > 1)]
  
  mat_t_df$name <- gsub(" ", ".", mat_t_df$name)
    
  mat_t_df <- merge(mat_t_df, dri, by.x = c("gene"), by.y = c("SYMBOL"))
  
  rownames(mat_t_df) <- mat_t_df$name
  
  dri_sel <- mat_t_df[, c("gene", "TIER")]
  
  mat_t_df <- mat_t_df[, !(names(mat_t_df) %in% c("chr_gene", "gene", "TIER", "count", "name"))]
  
  mat2plot <- as.matrix(t(mat_t_df), rownames.force = TRUE)
  
  
  ## construct heatmap
  
  options(device='x11', bitmapType="cairo")
    
  top_col <- list(Tier = c("High" = "darkseagreen2", "Mid" = "darkseagreen3", "Low" = "darkseagreen4"))
        
  top_ann <- HeatmapAnnotation(Tier = dri_sel$TIER[which(rownames(dri_sel) %in% colnames(mat2plot))],  col = top_col,
            annotation_legend_param = list(title = "Tier", at = c("High", "Mid", "Low"), labels = c("High", "Mid", "Low"),
                    title_gp = gpar(fontsize = 8), 
                    labels_gp = gpar(fontsize = 8))
  )
                    
  pal <- colorRampPalette(c("lightblue", "purple"))
  
  time_col <- pal(8)
            
  left_col <- list(Group = c("A1" = "gold2", "A2" = "darkorange2", "B" = "deeppink3"),
  
            Timepoint = c("PreTx" =  time_col[1], "ICI-4W" =  time_col[2], "C1D15" =  time_col[3], "C1W3" =  time_col[4], 
                    "C5D6" =  time_col[5], "PostTx" =  time_col[6], "LTSR" =  time_col[7], "PT" =  time_col[8]),
            
            Reference = c("buffycoat" = "darkred", 
                    "duodenum" = "deepskyblue", "duodenum-D2" = "deepskyblue3", "duodenum-D2-slow" = "deepskyblue4", 
                    "gastric" = "blue", "oesophagus" = "yellow", "PBMC" = "darkorchid"),
            
            Histology = c("EAC" = "plum1", "ESCC" = "darkblue"),
            
            Response = c("Responder" = "cyan", "Non-Responder" = "magenta")            
            
          )
  
  left_ann <- rowAnnotation(Group = metadata$Group, Timepoint = metadata$Time.Point, Reference = metadata$Reference, Histology = metadata$Histology, Response = metadata$Response, 
  
            col = left_col,
  
            annotation_legend_param = list(Group = list(title = "Group", at = c("A1", "A2", "B"), labels = c("A1", "A2", "B")),
            
              Timepoint = list(title = "Time Point", at = c("PreTx", "ICI-4W", "C1D15", "C1W3", "C5D6", "PostTx", "LTSR", "PT"), 
              
                labels = c("PreTx", "ICI-4W", "C1D15", "C1W3", "C5D6", "PostTx", "LTSR", "PT")),
                
              Reference = list(title = "Reference", at = c("buffycoat", "duodenum", "duodenum-D2", "duodenum-D2-slow", "gastric", "oesophagus", "PBMC"), 
              
                labels = c("buffycoat", "duodenum", "duodenum-D2", "duodenum-D2-slow", "gastric", "oesophagus", "PBMC")),
                
              Histology = list(title = "Histology", at = c("EAC", "ESCC"), labels = c("EAC", "ESCC")),
                
              Response = list(title = "Response", at = c("Responder", "Non-Responder"), labels = c("Responder", "Non-Responder"))              
              
            )
          )
            
            
  right_col <- list(QC = c("FAIL" = "red", "FLAG" = "darkorange", "PASS" = "green"),
  
            Gender = c("Male" = "forestgreen", "Female" = "blue"),
            
            Age = circlize::colorRamp2(c(20, max(as.numeric(metadata$Age), na.rm = TRUE)), c("lightblue", "purple"))
          )
  
  right_ann <- rowAnnotation(QC = metadata$QC, Gender = metadata$Gender, Age = metadata$Age, col = right_col,
  
            annotation_legend_param = list(QC = list(title = "QC", at = c("FAIL", "FLAG", "PASS"), labels = c("FAIL", "FLAG", "PASS")),
            
            Gender = list(title = "Gender", at = c("Female", "Male"), labels = c("Female", "Male")),
      
              Age = list(title = "Age", at = c(20, 40, 60, 80), labels = c("20", "40", "60", "80"))),
              
            show_legend = c(QC = TRUE, Gender = TRUE, Age = TRUE)
          )
            
  # from https://bioconductor.statistik.tu-dortmund.de/packages/3.1/bioc/vignettes/ComplexHeatmap/inst/doc/ComplexHeatmap.html
  # colors = structure(rainbow_hcl(4), names = c("1", "2", "3", "4"))
  
  # also see 
  # https://jokergoo.github.io/ComplexHeatmap-reference/book/more-examples.html
  # https://jokergoo.github.io/ComplexHeatmap-reference/book/legends.html#discrete-legends
  # https://www.nceas.ucsb.edu/sites/default/files/2020-04/colorPaletteCheatsheet.pdf
  
  column_labels = gsub(".", " ", colnames(mat2plot), fixed = TRUE)
    
  at <- c(0, 1, 2, 3, 4, 5)
  
  mat_col <- colorRamp2(c(0, 2.99, 3, 5), c("blue", "lightblue1", "wheat1", "red"))

  hmap <- Heatmap(matrix = mat2plot, 
        use_raster = T, raster_device = "png",
        col = mat_col,
        heatmap_legend_param = list(title = "CN", at = at, legend_gp = gpar(fill = mat_col(at)), color_bar = "discrete",
                        labels = c("homdel", "loh", "other loss", "no change", "gain", "big gain")),
        column_title = "Genes", row_title = "Samples",
        row_names_gp = gpar(fontsize = 7), # Text size for row names
        top_annotation = top_ann,
        right_annotation = right_ann,
        left_annotation = left_ann,
        column_labels = column_labels,
        column_names_gp = grid::gpar(fontsize = 8) # Text size for column names
      )
    
  pdf(file.path(segfile_dir, paste0(toTitleCase(segfile_name), "_", plot_name, ".pdf")), 
      width = 10, height = 6)
      
    draw(hmap, merge_legends = TRUE, heatmap_legend_side = "bottom")
    
  dev.off()
  

  png(file.path(segfile_dir, paste0(toTitleCase(segfile_name), "_", plot_name, "_mqc.png")), #type = "cairo", # maybe type = "cairo" or "cairo-png",
      #width = min(1440, (20 * ncol(mat2plot)) + 150), height = min(864, (10 * nrow(mat2plot)) + 150))
      width = 1440, height = 864)
      
    draw(hmap, merge_legends = TRUE, heatmap_legend_side = "bottom")
    
  dev.off()
    
  
  ## heatmap with alternative orientation
  
  mat_t <- t(mat2plot)
    
  bottom_col <- list(Group = c("A1" = "gold2", "A2" = "darkorange2", "B" = "deeppink3"),
  
            Timepoint = c("PreTx" =  time_col[1], "ICI-4W" =  time_col[2], "C1D15" =  time_col[3], "C1W3" =  time_col[4], 
                    "C5D6" =  time_col[5], "PostTx" =  time_col[6], "LTSR" =  time_col[7], "PT" =  time_col[8]),
            
            Reference = c("buffycoat" = "darkred", 
                    "duodenum" = "deepskyblue", "duodenum-D2" = "deepskyblue3", "duodenum-D2-slow" = "deepskyblue4", 
                    "gastric" = "blue", "oesophagus" = "yellow", "PBMC" = "darkorchid"),
            
            Histology = c("EAC" = "plum1", "ESCC" = "darkblue"),
            
            Response = c("Responder" = "cyan", "Non-Responder" = "magenta")
          )
  
  bottom_ann <- HeatmapAnnotation(Group = metadata$Group, Timepoint = metadata$Time.Point, Reference = metadata$Reference, Histology = metadata$Histology, Response = metadata$Response, col = bottom_col,
  
            annotation_legend_param = list(Group = list(title = "Group", at = c("A1", "A2", "B"), labels = c("A1", "A2", "B")),
            
              Timepoint = list(title = "Time Point", at = c("PreTx", "ICI-4W", "C1D15", "C1W3", "C5D6", "PostTx", "LTSR", "PT"), 
              
                labels = c("PreTx", "ICI-4W", "C1D15", "C1W3", "C5D6", "PostTx", "LTSR", "PT")),
                
              Reference = list(title = "Reference", at = c("buffycoat", "duodenum", "duodenum-D2", "duodenum-D2-slow", "gastric", "oesophagus", "PBMC"), 
              
                labels = c("buffycoat", "duodenum", "duodenum-D2", "duodenum-D2-slow", "gastric", "oesophagus", "PBMC")),
                
              Histology = list(title = "Histology", at = c("EAC", "ESCC"), labels = c("EAC", "ESCC")),
                
              Response = list(title = "Response", at = c("Responder", "Non-Responder"), labels = c("Responder", "Non-Responder"))
            )
          )
  
  # BA changed line below to accommodate alternative TIER values
  
  left_col <- list(Tier = c("High" = "darkseagreen2", "Mid" = "darkseagreen3", "Low" = "darkseagreen4"))
          
  left_ann <- rowAnnotation(Tier = dri_sel$TIER[which(rownames(dri_sel) %in% rownames(mat_t))],  col = left_col,
              annotation_legend_param = list(title = "Tier", at = c("High", "Mid", "Low"), labels = c("High", "Mid", "Low"),
                    title_gp = gpar(fontsize = 8), 
                    labels_gp = gpar(fontsize = 8))
          )
          
          
  top_col <- list(Gender = c("Male" = "forestgreen", "Female" = "blue"),
            Age = circlize::colorRamp2(c(20, max(as.numeric(metadata$Age), na.rm = TRUE)), c("lightblue", "purple")),
            
            QC = c("FAIL" = "red", "FLAG" = "darkorange", "PASS" = "green")
          )

  top_ann <- HeatmapAnnotation(Gender = metadata$Gender, Age = as.numeric(metadata$Age), QC = metadata$QC, col = top_col,
  
            annotation_legend_param = list(Gender = list(title = "Gender", at = c("Female", "Male"), labels = c("Female", "Male")),
      
              Age = list(title = "Age", at = c(20, 40, 60, 80), labels = c("20", "40", "60", "80")),
              
              QC = list(title = "QC", at = c("FAIL", "FLAG", "PASS"), labels = c("FAIL", "FLAG", "PASS"))),
              
            show_legend = c(Gender = TRUE, Age = TRUE, QC = TRUE)
          )
    
  mat_col <- colorRamp2(c(0, 2.99, 3, 5), c("blue", "lightblue1", "wheat1", "red"))
  
  at <- c(0, 1, 2, 3, 4, 5)
  
  hmap_t <- Heatmap(matrix = mat_t,
            use_raster = T, raster_device = "png",
            col = mat_col,
            heatmap_legend_param = list(title = "CN", at = 0:5, legend_gp = gpar(fill = mat_col(at)), color_bar = "discrete",
                        labels = c("homdel", "loh", "other loss", "no change", "gain", "big gain")),
            column_title = "Samples", row_title = "Genes",
            row_names_gp = gpar(fontsize = 7), # Text size for row names
            top_annotation = top_ann,
            left_annotation = left_ann,
            bottom_annotation = bottom_ann,
            row_labels = gsub(".", " ", rownames(mat_t), fixed = TRUE),
            column_names_gp = grid::gpar(fontsize = 8) # Text size for column names
          )
        
  pdf(file.path(segfile_dir, paste0(toTitleCase(segfile_name), "_transposed_", plot_name, ".pdf")), 
      width = 10, height = 6)
      
    draw(hmap_t, merge_legends = TRUE, heatmap_legend_side = "bottom")
    
  dev.off()
          
  png(file.path(segfile_dir, paste0(toTitleCase(segfile_name), "_", plot_name, "_transposed_mqc.png")), #type = "cairo", # maybe type = "cairo" or "cairo-png",
      #width = min(1440, (20 * ncol(mat_t)) + 150), height = min(864, (10 * nrow(mat_t)) + 150))
      width = 1440, height = 864)
      
    draw(hmap_t, merge_legends = TRUE, heatmap_legend_side = "bottom")
    
  dev.off()
  
  
  ## heatmap without dendrograms, ordered on sample pair and chr/gene for comparison between runs with different R_SEED
  
  top_col <- list(Tier = c("High" = "darkseagreen2", "Mid" = "darkseagreen3", "Low" = "darkseagreen4"))      
          
  top_ann <- HeatmapAnnotation(Tier = dri_sel$TIER[which(rownames(dri_sel) %in% colnames(mat2plot))],  col = top_col,
            annotation_legend_param = list(title = "Tier", at = c("High", "Mid", "Low"), labels = c("High", "Mid", "Low"),
                    title_gp = gpar(fontsize = 8), 
                    labels_gp = gpar(fontsize = 8))
  )
                    
  pal <- colorRampPalette(c("lightblue", "purple"))
  
  time_col <- pal(8)
            
  left_col <- list(Group = c("A1" = "gold2", "A2" = "darkorange2", "B" = "deeppink3"),
  
            Timepoint = c("PreTx" =  time_col[1], "ICI-4W" =  time_col[2], "C1D15" =  time_col[3], "C1W3" =  time_col[4], 
                    "C5D6" =  time_col[5], "PostTx" =  time_col[6], "LTSR" =  time_col[7], "PT" =  time_col[8]),
            
            Reference = c("buffycoat" = "darkred", 
                    "duodenum" = "deepskyblue", "duodenum-D2" = "deepskyblue3", "duodenum-D2-slow" = "deepskyblue4", 
                    "gastric" = "blue", "oesophagus" = "yellow", "PBMC" = "darkorchid"),
            
            Histology = c("EAC" = "plum1", "ESCC" = "darkblue"),
            
            Response = c("Responder" = "cyan", "Non-Responder" = "magenta")            
            
          )
  
  left_ann <- rowAnnotation(Group = metadata$Group, Timepoint = metadata$Time.Point, Reference = metadata$Reference, Histology = metadata$Histology, Response = metadata$Response, 
  
            col = left_col,
  
            annotation_legend_param = list(Group = list(title = "Group", at = c("A1", "A2", "B"), labels = c("A1", "A2", "B")),
            
              Timepoint = list(title = "Time Point", at = c("PreTx", "ICI-4W", "C1D15", "C1W3", "C5D6", "PostTx", "LTSR", "PT"), 
              
                labels = c("PreTx", "ICI-4W", "C1D15", "C1W3", "C5D6", "PostTx", "LTSR", "PT")),
                
              Reference = list(title = "Reference", at = c("buffycoat", "duodenum", "duodenum-D2", "duodenum-D2-slow", "gastric", "oesophagus", "PBMC"), 
              
                labels = c("buffycoat", "duodenum", "duodenum-D2", "duodenum-D2-slow", "gastric", "oesophagus", "PBMC")),
                
              Histology = list(title = "Histology", at = c("EAC", "ESCC"), labels = c("EAC", "ESCC")),
                
              Response = list(title = "Response", at = c("Responder", "Non-Responder"), labels = c("Responder", "Non-Responder"))              
              
            )
          )
            
            
  right_col <- list(QC = c("FAIL" = "red", "FLAG" = "darkorange", "PASS" = "green"),
  
            Gender = c("Male" = "forestgreen", "Female" = "blue"),
            
            Age = circlize::colorRamp2(c(20, max(as.numeric(metadata$Age), na.rm = TRUE)), c("lightblue", "purple"))
          )
  
  right_ann <- rowAnnotation(QC = metadata$QC, Gender = metadata$Gender, Age = metadata$Age, col = right_col,
  
            annotation_legend_param = list(QC = list(title = "QC", at = c("FAIL", "FLAG", "PASS"), labels = c("FAIL", "FLAG", "PASS")),
            
            Gender = list(title = "Gender", at = c("Female", "Male"), labels = c("Female", "Male")),
      
              Age = list(title = "Age", at = c(20, 40, 60, 80), labels = c("20", "40", "60", "80"))),
              
            show_legend = c(QC = TRUE, Gender = TRUE, Age = TRUE)
          )
            
  column_labels = gsub(".", " ", colnames(mat2plot), fixed = TRUE)
    
  at <- c(0, 1, 2, 3, 4, 5)
    
  mat_col <- colorRamp2(c(0, 2.99, 3, 5), c("blue", "lightblue1", "wheat1", "red"))
  
  hmap_ord <- Heatmap(matrix = mat2plot, 
        use_raster = T, raster_device = "png",
        col = mat_col,
        row_order = sort(rownames(mat2plot)),
        column_order = sort(colnames(mat2plot), method = "radix"),
        heatmap_legend_param = list(title = "CN", at = at, legend_gp = gpar(fill = mat_col(at)), color_bar = "discrete",
                        labels = c("homdel", "loh", "other loss", "no change", "gain", "big gain")),
        column_title = "Genes", row_title = "Samples",
        row_names_gp = gpar(fontsize = 7), # Text size for row names
        top_annotation = top_ann,
        right_annotation = right_ann,
        left_annotation = left_ann,
        column_labels = sort(gsub(".", " ", colnames(mat2plot), fixed = TRUE), method = "radix"),
        column_names_gp = grid::gpar(fontsize = 8) # Text size for column names  
      )
    
    pdf(file.path(segfile_dir, paste0(toTitleCase(segfile_name), "_", plot_name, "_ordered.pdf")), 
      width = 10, height = 6)
      
    draw(hmap_ord, merge_legends = TRUE, heatmap_legend_side = "bottom")
    
  dev.off()
  

  png(file.path(segfile_dir, paste0(toTitleCase(segfile_name), "_", plot_name, "_ordered_mqc.png")), #type = "cairo", # maybe type = "cairo" or "cairo-png",
      #width = min(1440, (20 * ncol(mat2plot)) + 150), height = min(864, (10 * nrow(mat2plot)) + 150))
      width = 1440, height = 864)
      
    draw(hmap_ord, merge_legends = TRUE, heatmap_legend_side = "bottom")
    
  dev.off()
}

