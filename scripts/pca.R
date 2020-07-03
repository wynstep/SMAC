################################################################################
#
#   File name: pca.R
#
#   Authors: Stefano Pirro', PhD ( s.pirro@qmul.ac.uk )
#
#   Barts Cancer Institute,
#   Queen Mary, University of London
#   Charterhouse Square, London EC1M 6BQ
#
#		Description: this script calculated the principal component analysis,
#		starting from an expression matrix
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
	source(paste0(GetScriptWD(), "/functions.R")) # file containing functions
	source(paste0(GetScriptWD(), "/vars.R")) # file containing vars
	# Taking advantage of the pacman module for installing unmet dependencies
	if (!require("pacman")) install.packages("pacman")
	pacman::p_load(dependencies[["pca"]], character.only = T)

#===============================================================================
#	Load Arguments for the analysis
#===============================================================================
	parser <- ArgumentParser()
	parser$add_argument("-e", "--exprFile", action="store", required=TRUE, dest="exprFile", help="Input expression file to analyse")
	parser$add_argument('-t', "--metaFile", action="store", required=TRUE, dest="metaFile", help="Metadata for the sample")
	parser$add_argument('-i', "--index", action="store", required=FALSE, type="integer", default=0, dest="index", help="File index")
	parser$add_argument("-n", "--nCores", action="store", required=TRUE, dest="nCores", type="integer", default=1, help="Number of available cores for the analysis")
	parser$add_argument("-g", "--nGenes", action="store", required=TRUE, dest="nGenes", type="integer", default=100, help="Number of genes to include in the analysis")
	parser$add_argument("-o", "--outputDir", action="store", required=TRUE, dest="outputDir", help="Output directory")
	args <- parser$parse_args()
	PC1 = 1
	PC2 = 2
	PC3 = 3
#===============================================================================
#	Load input files
#===============================================================================
	# load expression file
	exprData <- fread(args$exprFile, sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE, fill=TRUE, data.table=FALSE, nThread=args$nCores)
	# load metadata file
	metaData <- fread(args$metaFile, sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE, fill=TRUE, data.table=FALSE, nThread=args$nCores)
	# extract just values corrensponding to samples
		# Indexes of samples
	sIndex <- grep("GSM*", colnames(exprData))
		# Index of geneName
	gIndex <- grep("^gene(.+|)(symbol|name|id)", colnames(exprData))
		# Extract data
	exprData.numeric <- exprData[,sIndex]
		# set rownames (geneName)
	if (length(gIndex) > 0) {
		rownames(exprData.numeric) <- make.unique(exprData[,gIndex])
	} else {
		rownames(exprData.numeric) <- make.unique(as.character(exprData[,1]))
	}

		# Intersecting the samples between expression and target files
	selectedSamples <- intersect(metaData[,1],colnames(exprData.numeric))
		# Remove NAs rows and select samples present in the metaData
	exprData.numeric <- exprData.numeric[complete.cases(exprData.numeric),selectedSamples, drop=FALSE]

#===============================================================================
#	Map metadata info to file name
#===============================================================================
	# extract field informative
	toMatch <- c("source_name*","characteristics_*","title")
	for (tm in toMatch) {
		infField <- grep(tm, colnames(metaData))
		if (grep(tm, colnames(metaData)) > 0) { break }
	}
	metaSamples <- sapply(1:nrow(metaData), function(x)
									as.character(gsub("[][]", "",paste0(metaData[x,infField], collapse="--"))))
	metaSamples <- as.character(lapply(as.character(metaSamples), CutString, n_words=6)) # we cut the string if it's too long
	names(metaSamples) <- metaData[,1]

#===============================================================================
#     Performing PCA
#===============================================================================
	# Keep only probes with standard deviation > 0 across all samples
	rsd <- apply(exprData.numeric,1,sd)
	exprData.sd <- exprData.numeric[rsd>0,]

	# Perform principal components analysis
	exprData.pca <- prcomp(t(exprData.sd), scale=FALSE)

	# Get variance importance for all principal components
	importance.pca <- summary(exprData.pca)$importance[2,]
	# Here we transform the summary into a dataframe and we convert the variances into a percentage
	importance.pca.df <- data.frame(PC = names(importance.pca), Variances = as.numeric(round(summary(exprData.pca)$importance[2,]*100,2)), stringsAsFactors=FALSE)

#===============================================================================
#     Generate Plots
#===============================================================================

	#===============================================================================
	#     Barplot
	#===============================================================================

	# sorting data by value
	importance.pca.df$PC <- factor(importance.pca.df$PC, levels = importance.pca.df$PC)
	# Plotting
	p <- plot_ly(importance.pca.df, x = ~PC, y = ~Variances, type = 'bar', width = 800, height = 600) %>%
	layout(title = "The variances captured by principal components", xaxis = list(title = ""), margin = list(l=50, r=50, b=100, t=100, pad=4), autosize = F)
	# Saving as HTML
	widgetFile = paste(args$outputDir,paste("pcaBP",args$index,"html",sep="."),sep="/")
	saveWidgetFix(p, file=widgetFile)

	#===============================================================================
	#     Scatterplot
	#===============================================================================

	# Prepare DataFrame
	exprData.pca.df <- data.frame("Target" = metaSamples, "PC1" = exprData.pca$x[,PC1], "PC2" = exprData.pca$x[,PC2], "PC3" = exprData.pca$x[,PC3])

  # Plotting 2D scatterplot
  p <- plot_ly(exprData.pca.df, x = ~PC1, y = ~PC2, color = ~Target, text=rownames(exprData.pca.df), type='scatter', mode = "markers", marker = list(size=10, symbol="circle"), width = 800, height = 600) %>%
  layout(title = "", xaxis = list(title = paste("PC ", PC1, " (",importance.pca.df[PC1,"Variances"]," %)",sep="")), yaxis = list(title = paste("PC ", PC2, " (",importance.pca.df[PC2, "Variances"]," %)",sep="")),
		margin = list(l=50, r=50, b=50, t=50, pad=4), autosize = F, legend = list(orientation = 'h', y = 1.1))
  	# Saving as HTML
  widgetFile = paste(args$outputDir,paste("pca2D",args$index,"html",sep="."),sep="/")
  saveWidgetFix(p, file=widgetFile)

	# Plotting 3D scatterplot
	p <- plot_ly(exprData.pca.df, x = ~PC1, y = ~PC2, z = ~PC3, color = ~Target, text=rownames(exprData.pca.df), type='scatter3d', mode = "markers", marker = list(size=8, symbol="circle"), width = 800, height = 800) %>%
	layout(scene = list(xaxis = list(title = paste("PC ", PC1, " (",importance.pca.df[PC1, "Variances"]," %)",sep="")), yaxis = list(title = paste("PC ", PC2, " (",importance.pca.df[PC2, "Variances"]," %)",sep="")),
		zaxis = list(title = paste("PC ", PC3, " (",importance.pca.df[PC3, "Variances"]," %)",sep="")) ),
		margin = list(l=50, r=50, b=50, t=50, pad=4), autosize = F, legend = list(orientation = 'h', y = 1.1))
		# Saving as HTML
	widgetFile = paste(args$outputDir,paste("pca3D",args$index,"html",sep="."),sep="/")
	saveWidgetFix(p, file=widgetFile)

	# Close any open graphics devices
	graphics.off()
