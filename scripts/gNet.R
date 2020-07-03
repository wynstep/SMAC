################################################################################
#
#   File name: gNet.R
#
#   Authors: Stefano Pirro', PhD ( s.pirro@qmul.ac.uk )
#
#   Barts Cancer Institute,
#   Queen Mary, University of London
#   Charterhouse Square, London EC1M 6BQ
#
#		Description: this script creates a gene network (using mentha as reference)
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
	pacman::p_load(dependencies[["gNet"]], character.only = T)

#===============================================================================
#	Load Arguments for the analysis
#===============================================================================
	parser <- ArgumentParser()
	parser$add_argument("-i", "--geneFile", action="store", required=TRUE, dest="geneFile", help="Input gene list to analyse")
	parser$add_argument("-l", "--queryLabel", action="store", required=TRUE, dest="queryLabel", help="Label to assign to output file")
	parser$add_argument("-n", "--nCores", action="store", required=TRUE, dest="nCores", type="integer", default=1, help="Number of available cores for the analysis")
	parser$add_argument("-o", "--outputDir", action="store", required=TRUE, dest="outputDir", help="Output directory")
	args <- parser$parse_args()
		# load network file
	netFile = paste0(GetScriptWD(), "/../src/mentha.csv")

#===============================================================================
#	Load input files
#===============================================================================
		# load gene list
	geneData <- fread(args$geneFile, sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE, fill=TRUE, data.table=FALSE, nThread=args$nCores)
		# load network data
	netData <- fread(netFile, sep=";", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE, fill=TRUE, data.table=FALSE, nThread=args$nCores)
		# Filtering netData according to selected genes and thresholds (between 0.7 and 1)
	netDataFiltered <- netData[(netData[,2] %in% geneData[,2] | netData[,5] %in% geneData[,2]) & (netData[,7] >= 0.7 & netData[,7] <= 1), ]
		# Collecting total list of genes selected
	netGenes <- data.frame("geneName" = unique(c(netDataFiltered[,2], netDataFiltered[,5])), stringsAsFactors=FALSE)

		# Aggregate PMIDs and geneScores according to the geneName
	geneDataAggreg <- geneData %>%
			group_by(geneName) %>%
			mutate(PMID = paste(PMID, collapse=","), geneScore = paste(geneScore, collapse=",")) %>%
			distinct(geneName, PMID, geneScore)

		# Calculate the mean of geneScores from geneDataAggreg
	geneDataAggreg$meanScore <- sapply(strsplit((geneDataAggreg$geneScore), split=","), function(x) mean(as.numeric(x)))

#===============================================================================
#   Building network
#===============================================================================

	# mapping gene file with all the genes in the network
	mappedGenes = merge(netGenes, geneDataAggreg, by="geneName", all.x=T)
	# check if the merge succeed, otherwise treat genData as mapped
	if (nrow(mappedGenes) == 0) { mappedGenes = geneDataAggreg }

	#================================================
	# Creating Network data
	#================================================

	# creating genes (nodes) dataframe -- for visNetwork
	nodes <- data.frame(
		id=mappedGenes$geneName,
		exprValue = mappedGenes$meanScore,
		label=mappedGenes$geneName,
		color="black", # just a sample color
		group="gene",
		shape="circle",
		title = paste0("<a href='http://www.genecards.org/cgi-bin/carddisp.pl?gene=",mappedGenes$Gene,"' target=new><b>", mappedGenes$geneName ,"</b></a><br>
									 <b> Mean Enrichment score (PMIDs): </b>",round(mappedGenes$meanScore,digits=3),"<br>
									 <b> All Enrichment scores (PMIDs): </b>",mappedGenes$geneScore,"<br>
									 <b> PMIDs: </b>",
									lapply(strsplit(mappedGenes$PMID, split=","), function(x)
										paste0("<a href='https://pubmed.ncbi.nlm.nih.gov/",x,"' target=new><b>",x,"</b></a>",collapse=" "))),
		shadow = FALSE,
		stringsAsFactors = FALSE)

		# Coloring genes (nodes) according to the pubmed value
		nodesColoured <- ColourNodes(nodes)

		# creating edges -- for visNetwork
		edges <- data.frame(from = netDataFiltered[,2],
												to = netDataFiltered[,5],
												length = netDataFiltered[,7],
												title = paste0("<b> Source: </b>",netDataFiltered[,2]," <b> Target: </b>", netDataFiltered[,5] ,"<br>
																				<b> Interaction score: </b>",round(netDataFiltered[,7],digits=3),"")
												)
		#================================================
		# Saving Network
		#================================================
		widgetFile <- paste(args$outputDir,paste(args$queryLabel,"gNet","html",sep="."),sep="/")

		network <- visNetwork(nodesColoured, edges, height = "600px", width = "100%") %>%
			visEdges(smooth = TRUE, physics = TRUE) %>%
			visPhysics(stabilization = FALSE, solver = "barnesHut")
		# save newtwork into file
		saveWidgetFix(network, file = widgetFile, background = "white")
