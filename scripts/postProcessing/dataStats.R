################################################################################
#
#   File name: dataStats.R
#
#   Authors: Stefano Pirro', PhD ( s.pirro@mir-nat.com)
#
#
#		Description: this script calculate total numer of analysis in MAP and creates a plot
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
	source(paste0(GetScriptWD(), "/../functions.R")) # file containing functions
	source(paste0(GetScriptWD(), "/../vars.R")) # file containing vars
	# Taking advantage of the pacman module for installing unmet dependencies
	if (!require("pacman")) install.packages("pacman")
	pacman::p_load(dependencies[["dataStats"]], character.only = T)

#===============================================================================
#	Load Arguments for the analysis
#===============================================================================
	parser <- ArgumentParser()
	parser$add_argument("-i", "--mapFolder", action="store", required=TRUE, help="MAP data folder")
	parser$add_argument("-o", "--outputFile", action="store", required=TRUE, help="File where to store stats")
	args <- parser$parse_args()
  srcFolder <- paste(args$mapFolder,"src",sep="/")
  dataFolder <- paste(args$mapFolder,"data",sep="/")

#===============================================================================
#	Calculate number of papers
#===============================================================================
  # For calculating the number of papers, we count the unique PMIDs in the src folder
  literatureFile = paste0(srcFolder,"/literature.papersInfo.tsv")
  mirnasFile = paste0(srcFolder,"/mirnas.papersInfo.tsv")
  # load Files into datatables
  literatureData <- fread(literatureFile, header=TRUE, sep="\t", stringsAsFactors=F, fill=TRUE, data.table=FALSE)
  mirnasData <- fread(mirnasFile, header=TRUE, sep="\t", stringsAsFactors=F, fill=TRUE, data.table=FALSE)
  totalData <- rbind(literatureData,mirnasData)
  # We estimate the number of papers from Unique PMIDs
  nPapers <- length(unique(totalData$PMID))

#===============================================================================
#	Calculate number of analyses
#===============================================================================
  # Gene expression, gene networks and correlation are evaluated from the expression files
  nExprCmd <- paste0("find ",dataFolder," -type f -name \"*.exprs.*.tsv\" | wc -l")
  nExpr <- strtoi(system(nExprCmd, intern = TRUE))
  # PCA
  nPCACmd <- paste0("find ",dataFolder," -type f -name \"pca*.html\" | wc -l")
  nPCA <- strtoi(system(nPCACmd, intern = TRUE))/3
  # Molecular Classification
  nMclassCmd <- paste0("find ",dataFolder," -type f -name \"mClass.*.html\" | wc -l")
  nMClass <- strtoi(system(nMclassCmd, intern = TRUE))
  # Tumour Purity
  ntPurityCmd <- paste0("find ",dataFolder," -type f -name \"tPurity.*.html\" | wc -l")
  ntPurity <- strtoi(system(ntPurityCmd, intern = TRUE))
  # Receptor status
  nrStatusCmd <- paste0("find ",dataFolder," -type f -name \"rStatus.*.html\" | wc -l")
  nrStatus<- strtoi(system(nrStatusCmd, intern = TRUE))

#===============================================================================
#	Generate Dataframe
#===============================================================================
  dataStatsDF <- data.frame(
                    "Analysis" = c("Papers", "Gene expression", "Correlation", "Gene network", "PCA", "Molecular classification","Tumour purity","Receptor status"),
                    "Freq" = c(nPapers,nExpr,nExpr,nExpr,nPCA,nMClass,ntPurity,nrStatus)
                    )

#===============================================================================
#	Make plotly barplot
#===============================================================================

statsBP <- plot_ly(
              dataStatsDF,
              x = ~Analysis,
              y = ~Freq,
              color = ~Analysis,
              type = "bar"
            ) %>%
            layout(title = "", xaxis = list(title = ""), yaxis = list(title = "# Data", type = "log"), margin = list(l=50, r=50, b=100, t=50, pad=4), autosize = T, showlegend = FALSE)
saveWidgetFix(statsBP, file=args$outputFile)
