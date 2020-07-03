################################################################################
#
#   File name: updateTable.R
#
#   Authors: Stefano Pirro', PhD ( s.pirro@mir-nat.com)
#
#
#		Description: this script update the available analyses in a selected table, according to the presence
#   of key files
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
	pacman::p_load(dependencies[["updateTable"]], character.only = T)

#===============================================================================
#	Load Arguments for the analysis
#===============================================================================
	parser <- ArgumentParser()
	parser$add_argument("-a", "--papersFile", action="store", required=TRUE, help="TSV containing papers details")
	parser$add_argument("-b", "--genesFile", action="store", required=TRUE, help="TSV containing PMID/genes connections")
  parser$add_argument("-i", "--inputDir", action="store", required=TRUE, help="Directory where to check the presence of datasets")
	parser$add_argument("-o", "--outputDir", action="store", required=TRUE, help="Output directory")
	args <- parser$parse_args()

#===============================================================================
#	Load input data
#===============================================================================

  papersData <- fread(args$papersFile, sep="\t", header=TRUE, stringsAsFactors=FALSE, fill=TRUE, data.table=FALSE)
  genesData <- fread(args$genesFile, sep="\t", header=TRUE, stringsAsFactors=FALSE, fill=TRUE, data.table=FALSE)
  # Aggregate genesData dataframe according to the PMID column
  genesData.aggreg <- aggregate(genesData$geneName, list(genesData$PMID), paste, collapse=",")
  colnames(genesData.aggreg) <- c("PMID","geneName")
  # Merge PMID info and associated genes
  papersData.aggreg <- merge(papersData, genesData.aggreg, by="PMID", all.x=TRUE)

#===============================================================================
#	Recording available analyses
#===============================================================================

  # Adding the availableAnalyses column with all NAs (for now)
  papersData.aggreg$availableAnalyses <- "NA"
	# Replacing , with ; in selected columns
	papersData.aggreg$author <- gsub(',', ';', papersData.aggreg$author)
	papersData.aggreg$geneName <- gsub(',', ';', papersData.aggreg$geneName)
	papersData.aggreg$meshHeadings <- gsub(',', ';', papersData.aggreg$meshHeadings)

  # Iterate over each row in tableData and verify the presence of some keyFiles for identifying performed analyses
  for (i in 1:nrow(papersData.aggreg)) {
    # Logging operation
    cat("Parsing PMID ",i,"/",nrow(papersData.aggreg),"\n")
    # Initialise array that will contain the performed analysis
    perfAnalysis <- c()
    pmid <- papersData.aggreg[i,"PMID"]
    # Initialise PMID dir to look inside
    pmidDir <- paste(args$inputDir, pmid, sep="/")
    # We first check the presence of ExprData and metaData
    exprFiles <- list.files(pmidDir, pattern = glob2rx("*.exprs.*.tsv"), recursive=TRUE, full.names = TRUE)
    metaFiles <- list.files(pmidDir, pattern = glob2rx("*.meta.*.tsv"), recursive=TRUE, full.names = TRUE)
    pcaFiles <- list.files(pmidDir, pattern = glob2rx("pca*.html"), recursive=TRUE, full.names = TRUE)
    rStatusFiles <- list.files(pmidDir, pattern = glob2rx("rStatus.*.html"), recursive=TRUE, full.names = TRUE)
    tPurityFiles <- list.files(pmidDir, pattern = glob2rx("tPurity.*.html"), recursive=TRUE, full.names = TRUE)
    mClassFiles <- list.files(pmidDir, pattern = glob2rx("mClass.*.html"), recursive=TRUE, full.names = TRUE)
		# Enrichment files
		goFiles <- list.files(pmidDir, pattern = glob2rx("GeneOnthologyBP.*.png"), recursive=TRUE, full.names = TRUE)
		keggFiles <- list.files(pmidDir, pattern = glob2rx("KEGG.*.png"), recursive=TRUE, full.names = TRUE)
		reactomeFiles <- list.files(pmidDir, pattern = glob2rx("REACTOME.*.png"), recursive=TRUE, full.names = TRUE)
		if (length(goFiles) > 0 | length(keggFiles) > 0 | length(reactomeFiles) > 0) {
			print(goFiles);
			print(keggFiles);
			print(reactomeFiles);
		}

    if (length(exprFiles) > 0 & length(metaFiles) > 0) {
      perfAnalysis <- c(perfAnalysis, c("gene_expression","gene_network","correlation"))

      if (length(pcaFiles) > 0) { perfAnalysis <- c(perfAnalysis, "pca") }
      if (length(rStatusFiles) > 0) { perfAnalysis <- c(perfAnalysis, "receptor_status") }
      if (length(tPurityFiles) > 0) { perfAnalysis <- c(perfAnalysis, "tumour_purity") }
      if (length(tPurityFiles) > 0) { perfAnalysis <- c(perfAnalysis, "molecular_classification") }
			if (length(goFiles) > 0 | length(keggFiles) > 0 | length(reactomeFiles) > 0) { perfAnalysis <- c(perfAnalysis, "functional_enrichment") }

      # Update tableData with the value
      papersData.aggreg[i, "availableAnalyses"] <- paste(perfAnalysis, collapse=";")
    }
  }

#===============================================================================
#	Save new dataframe
#===============================================================================
  # Backup previous Table
system(paste("cp",args$papersFile,paste(args$papersFile,"bk",sep="."), sep=" "))
  # Create new table filename
newPapersFile <- paste(args$outputDir,basename(args$papersFile),sep="/")
  # Overwrite the table
fwrite(papersData.aggreg, file=newPapersFile, sep="\t", row.names=FALSE, col.names=TRUE, quote=TRUE)
