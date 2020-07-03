################################################################################
#
#   File name: gExpr.R
#
#   Authors: Stefano Pirro', PhD ( s.pirro@qmul.ac.uk )
#
#   Barts Cancer Institute,
#   Queen Mary, University of London
#   Charterhouse Square, London EC1M 6BQ
#
#		Description: this script plots an heatmap of the topN/-N differentially expressed genes
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
	pacman::p_load(dependencies[["gExpr"]], character.only = T)

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
	gIndex <- grep("^gene(.+|)(symbol|name|id)", colnames(exprData), ignore.case=TRUE)
		# Extract data
	exprData <- exprData[which(exprData[,gIndex] != ""), ]
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
	toMatch <- c("title", "characteristics_*")
	infField <- grep(paste(toMatch,collapse="|"), colnames(metaData))
	metaSamples <- sapply(1:nrow(metaData), function(x)
									as.character(gsub("[][]", "",paste0(metaData[x,infField], collapse="--"))))
	metaSamples <- as.character(lapply(as.character(metaSamples), CutString, n_words=6)) # we cut the string if it's too long
	names(metaSamples) <- metaData[,1]

#===============================================================================
#	Selection of differentially expressed genes
#===============================================================================
	# Calculate zscore of expression for 1 gene across the samples
	exprData.numeric.zscore <- as.data.frame(t(scale(t(exprData.numeric)))
	rownames(exprData.numeric.zscore) <- rownames(exprData.numeric)
	medianZscore <- as.numeric(unlist(apply(exprData.numeric.zscore, 1, function(x) median(x))))
	exprData.numeric.zscore$median <- medianZscore

	# sort according to the z-score of mean values and select top/bottom N genes
	exprData.numeric.zscore <- exprData.numeric.zscore[order(exprData.numeric.zscore$median, decreasing = TRUE),]
	# select top/bottom nGenes -- MEDIAN column is excluded
	nGenes = args$nGenes/2
	exprData.numeric.zscore.selected <- rbind(head(exprData.numeric.zscore[,-1], nGenes), tail(exprData.numeric.zscore[,-1], nGenes))
	# Get rid of too long gene names
	rownames(exprData.numeric.zscore.selected) <- make.unique(as.character(lapply(as.character(rownames(exprData.numeric.zscore.selected)), CutString, n_words=1)))

	# generates gene annotation
	topGenes.color <- c(rep("red",nGenes), rep("blue",nGenes))
	names(topGenes.color) <- c(rep("high",nGenes), rep("low",nGenes))

#===============================================================================
#     Expression heatmap
#===============================================================================
	# We take advantage of ComplexHeatmap for generating the plot

	# Assigning random color to any metaData category
	metaData.color <- getMetaColours(as.character(metaSamples), "discrete")

	# preparing topAnnotation (sample type)
	colAnn <- HeatmapAnnotation(samples=as.character(metaSamples),
														col=list(samples=metaData.color))
	# preparing rightAnnotation (up-down genes)
	rowAnn = rowAnnotation(df = data.frame(geneLevel=names(topGenes.color)),
													col=list(geneLevel = topGenes.color))

	# Preparing color for expression
	exprColor <- colorRamp2(c(min(exprData.numeric.zscore.selected), 0, max(exprData.numeric.zscore.selected)), c("blue", "white", "red"))

	rowText <- rowAnnotation(labels=anno_text(rownames(exprData.numeric.zscore.selected)),
													gp = gpar(col="black", cex=0.1, fontsize=10))

	#===============================================================================
	#     Heatmap
	#===============================================================================

	heatMap <- rowText + Heatmap(data.matrix(exprData.numeric.zscore.selected),
							cluster_rows = FALSE,
							show_row_names = FALSE,
							show_row_dend = FALSE,
							row_dend_reorder = FALSE,
							top_annotation = colAnn,
							row_names_side = "left",
							col = exprColor,
							row_names_gp = gpar(fontsize = 5, cex = 0.1),
				 			column_names_gp = gpar(fontsize = 15, cex = 0.3),
							heatmap_legend_param = list(title = "gEX (z-score)"),
							raster_device = "png",
							raster_quality = 10 ) + rowAnn

	# initialising file name
	heatMapFileName = paste(args$outputDir,paste("exprHM",args$index,"png",sep="."),sep="/")
	# initialising png file
	png(filename = heatMapFileName, bg = "white", width=1000, height=2000, units="px", res=150)
	# drawing heatmap
	draw(heatMap, annotation_legend_side = "bottom")
	# close device
	dev.off()

# Close any open graphics devices
graphics.off()
