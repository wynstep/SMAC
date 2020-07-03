################################################################################
#
#   File name: tPurity.R
#
#   Authors: Stefano Pirro', PhD ( s.pirro@qmul.ac.uk )
#
#   Barts Cancer Institute,
#   Queen Mary, University of London
#   Charterhouse Square, London EC1M 6BQ
#
#		Description: this script compute and plot the tumour purity for (cancer) samples
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
	pacman::p_load(dependencies[["tPurity"]], character.only = T)

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
	gIndex <- grep("^gene*", colnames(exprData))
		# Extract data
	exprData.numeric <- exprData[,sIndex]
		# set rownames (geneName)
	if (length(gIndex) > 0) {
		rownames(exprData.numeric) <- make.unique(exprData[,gIndex])
	} else {
		rownames(exprData.numeric) <- make.unique(as.character(exprData[,1]))
	}

#===============================================================================
#	Map metadata info to file name
#===============================================================================
	# extract field informative
	toMatch <- c("source_name*", "title", "characteristics_*")
	infField <- grep(paste(toMatch,collapse="|"), colnames(metaData))
	metaSamples <- sapply(1:nrow(metaData), function(x)
									as.character(gsub("[][]", "",paste0(metaData[x,infField], collapse="--"))))
	names(metaSamples) <- metaData[,1]

#===============================================================================
#	Selection of cancer samples
#===============================================================================
	# Selecting the samples where the number of "cancer" and "normal" words
	selMD <- countCancer(metaSamples, CancerToMatch, NormalToMatch)
	metaSamples.cancer <- as.character(lapply(as.character(metaSamples[selMD]), CutString, n_words=6)) # we cut the string if it's too long
	names(metaSamples.cancer) <- metaData[selMD,1]
	# Select expression data for cancer samples
	selectedSamples <- intersect(names(metaSamples.cancer),colnames(exprData.numeric))
	# Remove NAs columns and select samples present in the metaData
	exprData.numeric.cancer <- exprData.numeric[,selectedSamples, drop=FALSE]
	exprData.numeric.cancer <- exprData.numeric.cancer[,colSums(is.na(exprData.numeric.cancer))<nrow(exprData.numeric.cancer)]

#===============================================================================
#	Calculate tumour purity using the ESTIMATE package
#===============================================================================
	if (nrow(exprData.numeric.cancer) > 0 & ncol(exprData.numeric.cancer) > 0 & length(metaSamples.cancer) > 0) {
		# writing expression data filtered on a tmp file (needed by estimate package)
		tmpExprFileName <- paste(args$outputDir,paste("tmpExpr",args$index,"tsv",sep="."),sep="/")
		write.table(exprData.numeric.cancer,file=tmpExprFileName,row.names=T,col.names = T,sep="\t",quote=F)
		# filtering common genes using the built-in Estimate function
		estFile <- paste(args$outputDir,paste("CommonGenesEstimate",args$index,"gct",sep="."),sep="/")
		suppressMessages(filterCommonGenes(tmpExprFileName, estFile, id = "GeneSymbol"))
		# estimate purity score
		scoreFile <- paste(args$outputDir,paste("EstimateScore",args$index,"gct",sep="."),sep="/")
		suppressMessages(estimateScore(estFile, output.ds=scoreFile))

		#===============================================================================
		#     Scatterplot
		#===============================================================================

		# reading the file generated by Estimate and import into a dataframe
		estData=read.table(file=scoreFile,sep="\t",skip=2,header=T,row.names=1)
		estData.df = data.frame('File_name'=colnames(estData[,-1]),
														'Description'=sapply(colnames(estData[,-1]), function(x) metaSamples.cancer[[x]]),
														'StromalScore'=round(as.numeric(estData[1,-1]), digits = 4),
														'ImmuneScore'=round(as.numeric(estData[2,-1]), digits = 4),
														'ESTIMATEScore'=round(as.numeric(estData[3,-1]), digits = 4),
														'TumourPurity'=round(as.numeric(estData[4,-1]), digits = 4) )
		# Save into TSV file
		tPurityFileName <- paste(args$outputDir,paste("tPurity",args$index,"tsv",sep="."),sep="/")
		write.table(estData.df, file=tPurityFileName, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

		# Make 3D scatterplot
		p <- plot_ly( estData.df, x= ~StromalScore, y= ~ImmuneScore, z= ~ESTIMATEScore, color = ~TumourPurity,
									hoverinfo= 'text',
									text = paste('Sample Name:', estData.df$File_name,
															 '<br><br>StromalScore:', estData.df$StromalScore,
															 '<br>ImmuneScore:', estData.df$ImmuneScore,
															 '<br>ESTIMATEScore:', estData.df$ESTIMATEScore,
															 '<br><br><b>Tumour Purity</b>:', estData.df$TumourPurity),
									marker = EstimateCbStyle, width = '800px' ) %>%
									add_markers() %>%
									layout(scene = list(xaxis = list(title = 'StromalScore', orientation = "v"),
																		yaxis = list(title = 'ImmuneScore'),
																		zaxis = list(title = 'EstimateScore')),
									margin = plotMargins[["tPurity"]], title = "")

		# Saving widget
		widgetFile <- paste(args$outputDir,paste("tPurity",args$index,"html",sep="."),sep="/")
		saveWidgetFix(p, file=widgetFile)

		#===============================================================================
		#     Cleaning residues
		#===============================================================================

		# removing gct and score files
		file.remove(tmpExprFileName)
		file.remove(estFile)
		file.remove(scoreFile)

	}
