################################################################################
#
#   File name: PCA.R
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
suppressMessages(library(plotly))
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
  #make_option(c("-p", "--principal_component"), action="store", default=NA, type='character',
  #            help="The principal component to be plotted together with the two subsequent most prevalent principal components"),
  make_option(c("-d", "--dir"), action="store", default=NA, type='character',
              help="Default directory")
)

opt = parse_args(OptionParser(option_list=option_list))
expFile <- opt$exp_file
annFile <- opt$target
target <- opt$colouring
number <- opt$number
#PC1 <- as.numeric(opt$principal_component)
PC1 <- 1
PC2 <- PC1 + 1
PC3 <- PC1 + 2
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

  # subsetting expression file
  expData.subset <- expData[complete.cases(expData),colnames(expData) %in% selected_samples]

  # subsetting annData
  annData.subset <- annData[annData[,1] %in% selected_samples,]

  # getting index of target column in the annotation file
  target.index = grep(target, colnames(annData.subset))

  # creating mapping dictionary file_name -- group
  targets.groups = sapply(annData.subset[,1], function(x) {return(annData.subset[annData.subset[,1]==x,target.index])})

#===============================================================================
#     Principal components analysis
#===============================================================================

  # Keep only probes with standard deviation > 0 across all samples
  rsd <- apply(expData.subset,1,sd)
  expData.subset.sd <- expData.subset[rsd>0,]

  # Perform principal components analysis
  expData.pca <- prcomp(t(expData.subset), scale=FALSE)

  # Get variance importance for all principal components
  importance.pca <- summary(expData.pca)$importance[2,]
  # Here we transfomr the summary into a dataframe and we convert the variances into a percentage
  importance.pca.df <- data.frame(PC = names(importance.pca), Variances = as.numeric(round(summary(expData.pca)$importance[2,]*100,2)), stringsAsFactors=FALSE)

#===============================================================================
#     Generate Plots
#===============================================================================

  ### BarPlot ###
    # sorting data by value
    importance.pca.df$PC <- factor(importance.pca.df$PC, levels = importance.pca.df$PC)
    # Plotting
    p <- plot_ly(importance.pca.df, x = ~PC, y = ~Variances, type = 'bar', width = 800, height = 600) %>%
    layout(title = "The variances captured by principal components", xaxis = list(title = ""), margin = list(l=50, r=50, b=100, t=100, pad=4), autosize = F)
    # Saving as HTML
    widget_fn = paste(outFolder,paste0("pca_bp.",number,".html"),sep="/")
    suppressWarnings(saveWidgetFix(p, file=widget_fn, selfcontained = TRUE))

  ### ScatterPlots (3D, 2D) ###
    # Prepare DataFrame
	targets.labels <- as.character(lapply(as.character(targets.groups), CutString, n_words=2)) # we cut the string if it's too long
    expData.pca.df <- data.frame(targets.groups, expData.pca$x[,PC1], expData.pca$x[,PC2], expData.pca$x[,PC3])
    colnames(expData.pca.df) <- c("Target", "PC1", "PC2", "PC3")
    # Plotting 2D scatterplot
    p <- plot_ly(expData.pca.df, x = ~PC1, y = ~PC2, color = ~Target, text=rownames(expData.pca.df), type='scatter', mode = "markers", marker = list(size=10, symbol="circle"), width = 800, height = 600) %>%
    layout(title = "", xaxis = list(title = paste("PC ", PC1, " (",importance.pca.df[PC1,"Variances"]," %)",sep="")), yaxis = list(title = paste("PC ", PC2, " (",importance.pca.df[PC2, "Variances"]," %)",sep="")), margin = list(l=50, r=50, b=50, t=50, pad=4), autosize = F, legend = list(orientation = 'h', y = 1.1))
    # Saving as HTML
    widget_fn = paste(outFolder,paste0("pca_2d.",number,".html"),sep="/")
    suppressWarnings(saveWidgetFix(p, file=widget_fn, selfcontained = TRUE))
    # Plotting 3D scatterplot
    p <- plot_ly(expData.pca.df, x = ~PC1, y = ~PC2, z = ~PC3, color = ~Target, text=rownames(expData.pca.df), type='scatter3d', mode = "markers", marker = list(size=8, symbol="circle"), width = 800, height = 800) %>%
    layout(scene = list(xaxis = list(title = paste("PC ", PC1, " (",importance.pca.df[PC1, "Variances"]," %)",sep="")), yaxis = list(title = paste("PC ", PC2, " (",importance.pca.df[PC2, "Variances"]," %)",sep="")), zaxis = list(title = paste("PC ", PC3, " (",importance.pca.df[PC3, "Variances"]," %)",sep="")) ), margin = list(l=50, r=50, b=50, t=50, pad=4), autosize = F, legend = list(orientation = 'h', y = 1.1))
    # Saving as HTML
    widget_fn = paste(outFolder,paste0("pca_3d.",number,".html"),sep="/")
    suppressWarnings(saveWidgetFix(p, file=widget_fn, selfcontained = TRUE))

  # Close any open graphics devices
  graphics.off()
