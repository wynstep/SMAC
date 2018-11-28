################################################################################
#
#   File name: gene_expr.R
#
#   Authors: Stefano Pirro' ( s.pirro@qmul.ac.uk )
#
#   Barts Cancer Institute,
#   Queen Mary, University of London
#   Charterhouse Square, London EC1M 6BQ
#
################################################################################

#===============================================================================
#    Loading working directory
#===============================================================================
GetScriptWD <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(dirname(normalizePath(sub(needle, "", cmdArgs[match]))))
  } else {
    # 'source'd via R console
    return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  }
}

#===============================================================================
#    Load libraries
#===============================================================================
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))
suppressMessages(library(optparse))
suppressMessages(library(data.table))
source(paste0(GetScriptWD(), "/functions.R")) # file containing functions
source(paste0(GetScriptWD(), "/vars.R")) # file containing vars

#===============================================================================
#    Catching the arguments
#===============================================================================

option_list = list(
  make_option(c("-e", "--exp_file"), action="store", default=NA, type='character',
              help="File containing experimental data"),
  make_option(c("-t", "--target"), action="store", default=NA, type='character',
              help="Clinical data saved in tab-delimited format"),
  make_option(c("-c", "--colouring"), action="store", default=NA, type='character',
              help="Variable from the samples annotation file to be used for samples colouring"),
  make_option(c("-n", "--number"), action="store", default=NA, type='character',
              help="Number of the plot"),
  make_option(c("-d", "--dir"), action="store", default=NA, type='character',
              help="Default directory")
)

opt = parse_args(OptionParser(option_list=option_list))
expFile <- opt$exp_file
annFile <- opt$target
target <- opt$colouring
number <- opt$number
outFolder <- opt$dir


#===============================================================================
#   Pre-processing of the input files
#===============================================================================

  # Read file with expression data
  expData <- read.table(expFile,sep="\t",header=TRUE, stringsAsFactors = FALSE, row.names=NULL)
  rownames(expData) = make.unique(expData[,1])

  # Read sample annotation file
  annData <- fread(annFile,sep="\t",header=TRUE, data.table=FALSE)

  # Intersecting the samples between expression and target files
  selected_samples <- intersect(as.character(annData[,1]),colnames(expData))

  # subsetting expression file by selected samples
  expData.subset <- expData[,colnames(expData) %in% selected_samples]

  # subsetting annData
  annData.subset <- annData[annData[,1] %in% selected_samples,]

  # getting index of target column in the annotation file
  target.index = grep(target, colnames(annData.subset))

  # getting groups names
  groups = unique(annData.subset[,target.index])

  # creating mapping dictionary file_name -- group
  targets.groups = sapply(annData.subset[,1], function(x) {return(annData.subset[annData.subset[,1]==x,target.index])})

#===============================================================================
#     Expression heatmap (top 20, 50, 100 up-regulated/down-regulated genes across all samples)
#===============================================================================
  
  # transform expression matrix into a z-score matrix (by row)
  expData.scaled <- as.data.frame(t(scale(t(expData.subset))))
  # extract mean z-score for each row (gene), then desc sort
  genes.zscore <- sort(rowMeans(expData.scaled), decreasing = T)
  
  # Getting the top 20, 50, 100 top/down regulated genes (across all samples)
  for (i in c(20, 50, 100)){
    # selecting the top highest expressed
    top_genes = genes.zscore[1:i]
    # assigning red color to the category
    top_genes_up_down = rep("red", i)
    # selecting the lowest expressed
    top_genes = c(top_genes, tail(genes.zscore, i))
    # assigning blue color to the category
    top_genes_up_down = c(top_genes_up_down, rep("blue", i))
    names(top_genes_up_down) = c(rep("high", length(top_genes)/2), rep("low", length(top_genes)/2))
    # subset the expression matrix (normalised)
    expData.scaled.uselected = expData.scaled[names(top_genes),]
    # color target groups by value
    targets.groups.color = getTargetsColours(as.character(targets.groups))
    # Plotting the heatmap with extracted values
      # top annotation (sample type)
      col_ann = HeatmapAnnotation(df = data.frame(sample_type=targets.groups), 
                                  col=list(sample_type = targets.groups.color)
                                  )
      # right annotation (up - down genes)
      row_ann = rowAnnotation(df = data.frame(gene_level=names(top_genes_up_down)), 
                              col=list(gene_level = top_genes_up_down)
                              )
      row_text = rowAnnotation(labels=anno_text(names(top_genes), 
                               which = "row",
                               gp = gpar(col="black", cex=0.5)))
      
      # Plotting
       hm <- row_text + Heatmap(expData.scaled.uselected,
                     cluster_rows = FALSE,
                     show_row_names = FALSE,
                     show_row_dend = FALSE,
                     row_dend_reorder = FALSE,
                     top_annotation = top_ann,
                     row_names_side = "left",
                     heatmap_legend_param = list(title = "Expression scaled (z-score)"),
                     raster_device = "png",
                     raster_quality = 10
                    ) + row_ann 
       
       # initialising file name
       hm_fn = paste(outFolder,paste0("expression_heatmap_scaled_values.",i,".png"),sep="/")
       # initialising png file
       png(filename = hm_fn, width = 4000, height = 30*nrow(expData.scaled.uselected), units = "px", bg = "white", res = 300)
       # drawing heatmap
       draw(hm)
       # close device
       dev.off()
  }
      