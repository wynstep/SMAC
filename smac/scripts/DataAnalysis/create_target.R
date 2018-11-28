################################################################################
#
#   File name: create_target.R
#   Description: R script to create a target file for maximising the information rate
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
suppressMessages(library(heatmaply))
suppressMessages(library(optparse))
suppressMessages(library(data.table))
source(paste0(GetScriptWD(), "/functions.R")) # file containing functions
source(paste0(GetScriptWD(), "/vars.R")) # file containing vars

#===============================================================================
#    Catching the arguments
#===============================================================================
option_list = list(
  make_option(c("-p", "--pheno"), action="store", default=NA, type='character',
              help="Clinical data saved in tab-delimited format"),
  make_option(c("-c", "--cont"), action="store", default=NA, type='character',
              help="Index of the phenotypic data"),
  make_option(c("-d", "--dir"), action="store", default=NA, type='character',
              help="Default directory")
)

opt = parse_args(OptionParser(option_list=option_list))
pFile <- opt$pheno
cont <- opt$cont
outFolder <- opt$dir

#===============================================================================
#   Pre-processing of the input files
#===============================================================================

# Read sample annotation file
pData <- fread(pFile,sep="\t",header=TRUE, data.table=FALSE)

# Get index of sample names
sample_name_index <- grep("geo.+accession", names(pData))

#===============================================================================
#     Compute "information score" for the phenoData
#===============================================================================

# Generate dataframe with the scores
pData.score <- apply(pData, c(1, 2), function(x) computeInfoScore(x))
# Calculare the avergae for each column
pData.average.score <- as.data.frame(colMeans(pData.score))
# Retrieve column with maximum value
max.score.column = which(pData.average.score==max(pData.average.score))[1] # in case of multiple max columns, take the first column in the order

#===============================================================================
#     Generate Plots
#===============================================================================

  ### Heatmap ###
    # Plotting
    pData.score[pData.score == 0] <- NA
    p <- heatmaply(pData.score, dendrogram="none", colors = colorRampPalette(c('dark blue','white','dark red'))(100), scale="none", main = "Info score", trace="none", hide_colorbar = FALSE, fontsize_row = 8, fontsize_col = 8) %>%
    layout(autosize = TRUE, width = 800, margin = list(l=150, r=50, b=150, t=50, pad=4), showlegend = FALSE)

    # Saving as HTML
    widget_fn = paste(outFolder,paste0("info_score.",cont,".html"),sep="/")
    saveWidgetFix(p, file=widget_fn, selfcontained = TRUE)

#===============================================================================
#     Save target file
#===============================================================================

  target <- data.frame("sample_name" = pData[,sample_name_index], "target" = pData[,max.score.column])
  OutFile = paste(outFolder,paste0("target.",cont,".tsv"),sep="/")
  write.table(target, file=OutFile, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
