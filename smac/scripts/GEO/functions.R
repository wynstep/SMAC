################################################################################
#
#   File name: functions.R
#   Description: This file contains all the important functions for GEO section
#
#   Authors: Stefano Pirro' ( s.pirro@qmul.ac.uk )
#
#   Barts Cancer Institute,
#   Queen Mary, University of London
#   Charterhouse Square, London EC1M 6BQ
#
################################################################################

#===============================================================================
#    Load libraries
#===============================================================================
suppressMessages(library(data.table))

# This function reduce gene redundancy by selecting only the most variable genes
selectDuplicatedGenes <- function(expData) {
  # calculate IQR for each gene in the expression file
  expData$var <- apply(expData[,2:ncol(expData)],1,function(x) IQR(x,na.rm=FALSE, type=7))
  # slimming expData by selecting the most variable genes
  expData.slimmed <- as.data.frame(setDT(expData)[, .SD[which.max(var)], by=eval(colnames(expData)[1])])
  # Making adjustments to the expData.slimmed
  rownames(expData.slimmed) = expData.slimmed[,1] # setting rownames
  expData.slimmed[,c(1,ncol(expData.slimmed))] <- NULL # removing first column (with GeneName) and last column (var)

  return(expData.slimmed)
}
