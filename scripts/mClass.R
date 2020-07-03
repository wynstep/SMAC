################################################################################
#
#   File name: mClass.R
#
#   Authors: Stefano Pirro', PhD ( s.pirro@qmul.ac.uk )
#
#   Barts Cancer Institute,
#   Queen Mary, University of London
#   Charterhouse Square, London EC1M 6BQ
#
#		Description: this script compute and plot the molecular classification for (cancer) samples
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
	pacman::p_load(dependencies[["mClass"]], character.only = T)

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
	# load pam50 genes
	pam50File = paste0(GetScriptWD(), "/../src/pam50.txt")
	pam50Data <- fread(pam50File, sep="\t", header=TRUE, stringsAsFactors= FALSE, check.names=FALSE, fill=TRUE, data.table=FALSE, nThread=args$nCores)
	# load expression file
	exprData <- fread(args$exprFile, sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE, fill=TRUE, data.table=FALSE, nThread=args$nCores)
	# load metadata file
	metaData <- fread(args$metaFile, sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE, fill=TRUE, data.table=FALSE, nThread=args$nCores)
	# extract just values corrensponding to samples
		# Indexes of samples
	sIndex <- grep("GSM*", colnames(exprData))
		# Index of geneName
	gIndex <- grep("^gene*", colnames(exprData), ignore.case=TRUE)
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
	# Remove NAs rows and select samples present in the metaData
	exprData.numeric.cancer <- exprData.numeric[,selectedSamples, drop=FALSE]
	exprData.numeric.cancer <- exprData.numeric.cancer[,colSums(is.na(exprData.numeric.cancer))<nrow(exprData.numeric.cancer)]

#===============================================================================
#	Detect pam50 subtypes
#===============================================================================
	if (nrow(exprData.numeric.cancer) > 0 & ncol(exprData.numeric.cancer) > 0 & length(metaSamples.cancer) > 0) {
		# Performing class prediction
		mClassPred <- calculateMclass(exprData.numeric.cancer, pam50Data)
		mClassReport <- data.frame("SampleName" = names(mClassPred$subtype),
															"Subtype" = mClassPred$subtype,
															"Description" = sapply(names(mClassPred$subtype), function(x) metaSamples[[x]]))
		# Save into TSV file
		mClassFileName <- paste(args$outputDir,paste("mClass",args$index,"tsv",sep="."),sep="/")
		write.table(mClassReport, file=mClassFileName, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


		# We generate the plot just if mClassReport is not empty
		if (nrow(mClassReport) > 0) {

			#===============================================================================
			#     Barplot
			#===============================================================================

			# calculating frequency of different molecular subtypes
			mClassFreq <- as.data.frame(table(mClassReport$Subtype))
			# calculating the total number of samples into the report
			totalSamples = sum(mClassReport$Freq)

			# Plotting
			p <- plot_ly(mClassFreq, x = ~Var1, y = ~Freq, type = 'bar',
							marker = list(color = c('rgba(255,51,51,0.3)', 'rgba(0,0,204,0.3)',
																			'rgba(102,204,0,0.3)', 'rgba(0,102,0,0.3)',
																			'rgba(204,0,102,0.3)'),
							# defining bar line colors (same as bar colors but different opacity)
							line = list(color = c('rgba(255,51,51,1)', 'rgba(0,0,204,1)',
																		'rgba(102,204,0,1)', 'rgba(0,102,0,1)',
																		'rgba(204,0,102,1)'),
							width = 1.5) ),
						# adding percentage in the hover label text
						text = paste(round((mClassFreq$Freq/totalSamples)*100, digits = 2),'%'),
						textpositionsrc = 'center',
						width = '800px' ) %>%
					layout(yaxis = list(title = 'Frequency'), xaxis = list(title = 'Subtypes'), barmode = 'group',
								title = 'Molecular classification', showlegend = FALSE)

			# Saving as HTML
			widgetFile <- paste(args$outputDir,paste("mClass",args$index,"html",sep="."),sep="/")
			saveWidgetFix(p, file=widgetFile)

		}
	}
