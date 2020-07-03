################################################################################
#
#   File name: rra.R
#
#   Authors: Stefano Pirro', PhD ( s.pirro@qmul.ac.uk )
#
#   Barts Cancer Institute,
#   Queen Mary, University of London
#   Charterhouse Square, London EC1M 6BQ
#
#		Description: this script performs a Robust Ranking Aggregation
#		on MeSH terms, sorted in a different way
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
	pacman::p_load(dependencies[["rra"]], character.only = T)
	# Load custom version of RRA package
	#for (f in list.files(path=paste0(GetScriptWD(), "/../src/R/RankAggreg/"), pattern="*.R", full.names = TRUE)) { source(f) }

#===============================================================================
#	Load Arguments for the analysis
#===============================================================================
	parser <- ArgumentParser()
	parser$add_argument("-a", "--meshLevel", action="store", required=TRUE, dest="meshLevel", help="List of MeSH sorted by level")
	parser$add_argument("-b", "--meshGAboundance", action="store", required=TRUE, dest="meshGAboundance", help="List of MeSH sorted by global aboundance")
	parser$add_argument("-c", "--meshLAboundance", action="store", required=TRUE, dest="meshLAboundance", help="List of MeSH sorted by local aboundance")
	parser$add_argument("-l", "--queryLabel", action="store", required=TRUE, dest="queryLabel", help="Label to assign to output file")
	parser$add_argument("-n", "--nCores", action="store", required=TRUE, dest="nCores", type="integer", default=1, help="Number of available cores for the analysis")
	parser$add_argument("-g", "--graphDir", action="store", required=TRUE, dest="graphDir", help="Directory to store intermediate graphs")
	parser$add_argument("-o", "--outputDir", action="store", required=TRUE, dest="outputDir", help="Output directory")
	args <- parser$parse_args()

#===============================================================================
#	Load input files
#===============================================================================
	# load meshLevel file
	meshLevelData <- fread(args$meshLevel, sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE, fill=TRUE, data.table=FALSE, nThread=args$nCores, showProgress=TRUE)
	# load meshGAboundance file
	meshGAboundanceData <- fread(args$meshGAboundance, sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE, fill=TRUE, data.table=FALSE, nThread=args$nCores, showProgress=TRUE)
	# load meshLAboundance file
	meshLAboundanceData <- fread(args$meshLAboundance, sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE, fill=TRUE, data.table=FALSE, nThread=args$nCores, showProgress=TRUE)

#===============================================================================
#	Find common mesh terms and subset MeSH data
#===============================================================================

	meshData.merged <- Reduce( function(x, y) merge(x, y, by="mesh"), list(meshLevelData, meshGAboundanceData, meshLAboundanceData))
	meshData.merged$mesh <- make.unique(meshData.merged$mesh) # Avoiding the presence of duplicates

#===============================================================================
#	Compute local/global aboundance
#===============================================================================
	meshData.merged$abRatio <- round(meshData.merged$lAboundance/meshData.merged$gAboundance, digits=2)

#===============================================================================
#	Sort MeSH terms (separately)
#===============================================================================
	meshLevelData.sorted = as.character(meshData.merged[order(-meshData.merged$depthLevel),"mesh"]) # descending
	meshLAboundanceData.sorted = as.character(meshData.merged[order(-meshData.merged$lAboundance),"mesh"]) # descending
	meshGAboundanceData.sorted = as.character(meshData.merged[order(meshData.merged$gAboundance),"mesh"]) # ascending
	meshRatio.sorted = as.character(meshData.merged[order(-meshData.merged$abRatio),"mesh"]) # descending

#===============================================================================
#	Aggregate sorted lists
#===============================================================================
	mergedMeSH <- list(meshLevelData.sorted,
										meshLAboundanceData.sorted,
										meshGAboundanceData.sorted,
										meshRatio.sorted)

	# Apply Rank Aggregation Method
	AggregatedList <- aggregateRanks(glist = mergedMeSH, full = TRUE, exact = TRUE)
	# Save into dataframe
	AggregatedListFileName <- paste(args$outputDir,paste(args$queryLabel,"AggregatedMeSH","tsv",sep="."),sep="/")
	write.table(AggregatedList, file = AggregatedListFileName, sep = "\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

#===============================================================================
#	Calculate p-value /score fluctuations
#===============================================================================

	# creating ordered list of positions (first vector)
query <- c(1:nrow(AggregatedList))
	# iterating 1000 times -- parameter can be changed later
	# create tau and p-value dataframes
tauContainer = data.frame(iteration=numeric(), tau=numeric())
pvalueContainer = data.frame(iteration=numeric(), pvalue=numeric())

for (i in 1:1000) {
  #Sys.sleep(0.1)
  cat("Calculating tau and p-value for iteration ", i, "\n")
  # create vector of random positions (from mesh dictionary, as background)
  background <- sample(x = c(1:allMeshLength), size = length(query), replace = FALSE)
  # calculate the kendall distance and the p-value
  distance <- cor.test(query,background,method = "kendall")
  # report and store tau distance and p-values (for each iteration)
  tauContainer[i,] <- c(i, distance$estimate)
  pvalueContainer[i,] <- c(i, distance$p.value)
}

# create plots for p-value, tau and rho fluctuations

# create plots for p-value and tau fluctuations
png(filename = paste(args$graphDir,paste(args$queryLabel,"tauFluctuations","png",sep="."),sep="/"))
  meanTau = round(mean(tauContainer$tau), digits = 3)
  tauPlot <- ggplot(data=tauContainer, aes(x=iteration, y=tau, group=1, colour = tau)) + geom_line(alpha = 0.1) + geom_point(alpha = 0.2)
  tauPlot + scale_color_gradient(low="blue", high="red") +
  geom_hline(yintercept = meanTau) +
  annotate("text", min(tauContainer$iteration+200), meanTau, vjust = -1, label = meanTau)
dev.off()

png(filename = paste(args$graphDir,paste(args$queryLabel,"pvalueFluctuations","png",sep="."),sep="/"))
  meanPvalue = round(mean(pvalueContainer$pvalue), digits = 3)
  pvaluePlot <- ggplot(data=pvalueContainer, aes(x=iteration, y=pvalue, group=1, colour = pvalue)) + geom_line(alpha = 0.1) + geom_point(alpha = 0.2)
  pvaluePlot + scale_color_gradient(low="blue", high="red") +
  geom_hline(yintercept = meanPvalue) +
  annotate("text", min(pvalueContainer$iteration+200), meanPvalue, vjust = -1, label = meanPvalue)
dev.off()

png(filename = paste(args$graphDir,paste(args$queryLabel,"rhoFluctuations","png",sep="."),sep="/"))
  meanRho = round(mean(AggregatedList$Score), digits = 3)
  rhoPlot <- ggplot(data=AggregatedList, aes(x=c(1:nrow(AggregatedList)), y=Score, group=1, colour = Score)) + geom_line(alpha = 0.1) + geom_point(alpha = 0.2)
  rhoPlot + scale_color_gradient(low="blue", high="red") +
  geom_hline(yintercept = meanRho) +
  annotate("text", min(nrow(AggregatedList)+200), meanRho, vjust = -1, label = meanRho)
dev.off()

# save p_values and fluctuations details
write.table(pvalueContainer, file = paste(args$outputDir,paste(args$queryLabel,"pvalueFluctuations","csv",sep="."),sep="/"), sep = ",", quote = FALSE, col.names = T, row.names = F)
write.table(tauContainer, file = paste(args$outputDir,paste(args$queryLabel,"tauFluctuations","csv",sep="."),sep="/"), sep = ",", quote = FALSE, col.names = T, row.names = F)
