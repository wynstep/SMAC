################################################################################
#
#   File name: rStatus.R
#
#   Authors: Stefano Pirro', PhD ( s.pirro@qmul.ac.uk )
#
#   Barts Cancer Institute,
#   Queen Mary, University of London
#   Charterhouse Square, London EC1M 6BQ
#
#		Description: this script compute and plot the receptor status for (cancer) samples
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
	pacman::p_load(dependencies[["rStatus"]], character.only = T)

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

	#===============================================================================
	#	This is a very important step -- Selecting the receptor genes from the dataframe
	#===============================================================================
	exprData.numeric.cancer <- exprData.numeric[rownames(exprData.numeric) %in% receptorGenes, selectedSamples, drop=FALSE]
	# Remove samples where all the 3 receptor values are NA
	exprData.numeric.cancer <- exprData.numeric.cancer[,colSums(is.na(exprData.numeric.cancer))<nrow(exprData.numeric.cancer)]
	# Refine selectedSamples
	selectedSamples <- intersect(selectedSamples, colnames(exprData.numeric.cancer))
	# If the receptor genes are not present in the dataframe, we simply cannot continue
	if (nrow(exprData.numeric.cancer) == length(receptorGenes) & ncol(exprData.numeric.cancer) > 0 & length(metaSamples.cancer) > 0) {

		#===============================================================================
		#	Calculate receptor status for cancer samples
		#===============================================================================
				# creating receptor status report dataframe
				rstatus.report <- calculateRstatus(exprData.numeric.cancer)
				rstatus.report.file <- as.data.frame(t(rstatus.report))
				# add sample names and description to final dataframe report
				rstatus.report.file$FileName = rownames(rstatus.report.file)
				rstatus.report.file$Description = sapply(rownames(rstatus.report.file), function(x) metaSamples.cancer[[x]])
				# write final rstatus report
				rStatusFileName = paste(args$outputDir,paste("rStatus",args$index,"tsv",sep="."),sep="/")
				write.table(rstatus.report.file, rStatusFileName, sep="\t", row.names = FALSE, col.names=TRUE, quote = FALSE)

			#===============================================================================
			#     BarPlot
			#===============================================================================

				# Important note! The triple negative status can be assessed JUST if all the status
				# for ALL the 3 markers (ER, PR, HER2), is present

				# Preparing dataframe for plotting
				rstatus.report.freq <- as.data.frame(apply(rstatus.report, 1, function(x) table(x)))
				# Preparing widget file
				widgetFile <- paste(args$outputDir,paste("rStatus",args$index,"html",sep="."),sep="/")

				# In all cases, we plot the status for individual markers
				p.status <- plot_ly(rstatus.report.freq,
												x = colnames(rstatus.report.freq),
												y =  as.numeric(rstatus.report.freq[1,]), # first row contains negative values
												type = 'bar',
												name = 'Negative',
												legendgroup="1",
												# rgba(255,51,51) --> red
												marker = list(color = 'rgba(255,51,51,0.3)',
																			line = list(color = 'rgba(255,51,51,1)',
																									width = 1.5)),
												text = paste(round((rstatus.report.freq[1,]/sum(rstatus.report.freq[,1]))*100, digits = 2),'%'),
												textpositionsrc = 'center',
												width = '800px' ) %>%
										add_trace(rstatus.report.freq,
												x = colnames(rstatus.report.freq),
												y = as.numeric(rstatus.report.freq[2,]), # second row contains positive values
												type = 'bar',
												name = 'Positive',
												legendgroup="1",
												# rgba(0,0,204) --> blue
												marker = list(color = 'rgba(0,0,204,0.3)',
																			line = list(color = 'rgba(0,0,204,1)',
																									width = 1.5)),
												text = paste(round((rstatus.report.freq[2,]/sum(rstatus.report.freq[,1]))*100, digits = 2),'%'),
												textpositionsrc = 'center',
												width = '800px' ) %>%
										layout(yaxis = list(title = 'Count'),
													 barmode = 'stack')

				#===============================================================================
		 		#     Plot Triple Negative
		 		#===============================================================================

				if (nrow(rstatus.report) == 3) { # this means we have the status for all the markers...
					# Here we count the number of triple negative samples
					tn.samples = sum(colSums(rstatus.report) == 0)
					ntn.samples = sum(colSums(rstatus.report) > 0)
					tot = tn.samples + ntn.samples

					p.tneg <- plot_ly(rstatus.report.freq,
														x = "TripleNegative",
														y = ntn.samples ,
														type = 'bar',
														name = 'No',
														legendgroup="2",
														marker = list(color = 'rgba(102,204,0,0.3)', # light green
																					line = list(color = 'rgba(102,204,0,1)', width = 1.5)),
																					text = paste(round((ntn.samples/tot)*100, digits = 2),'%'),
																					textpositionsrc = 'center',
																					width = '800px' ) %>%
													add_trace(rstatus.report.freq,
														x = "TripleNegative",
														y = tn.samples,
														type = 'bar',
														name = 'Yes',
														legendgroup="2",
														marker = list(color = 'rgba(204,0,102,0.3)', # dark pink
																					line = list(color = 'rgba(204,0,102,1)', width = 1.5)),
						 															text = paste(round((tn.samples/tot)*100, digits = 2),'%'),
						 															textpositionsrc = 'center',
						 															width = '800px' ) %>%
													layout(yaxis = list(title = 'Count'), barmode = 'stack')

													# merging plots together
													rstatus.plot <- subplot(p.status, p.tneg, nrows = 1, margin = 0.02, shareY = TRUE)
													rstatus.plot <- layout(rstatus.plot, title="Receptor Status");

					# Save widget
					saveWidgetFix(rstatus.plot, file=widgetFile)
				} else {
					# Save widget (without TN status)
					saveWidgetFix(p.status, file=widgetFile)
				}

	} else {
		print("Sorry, no receptor genes in the dataframe, cannot calculate receptor status")
	}
