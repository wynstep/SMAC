################################################################################
#
#   File name: DownloadGSE.R
#	  Description: This script connects to GEO and downloads all the relevant data
#
#   Authors: Stefano Pirro' ( s.pirro@qmul.ac.uk )
#
#   Barts Cancer Institute,
#   Queen Mary, University of London
#   Charterhouse Square, London EC1M 6BQ
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
suppressMessages(library(optparse))
suppressMessages(library(GEOquery))
source(paste0(GetScriptWD(), "/functions.R")) # file containing functions

#===============================================================================
#    Catching the arguments
#===============================================================================

option_list = list(
  make_option(c("-g", "--gse"), action="store", default=NA, type='character',
              help="GSE dataset to download"),
  make_option(c("-d", "--dir"), action="store", default=NA, type='character',
              help="Directory where to save data")
)

opt = parse_args(OptionParser(option_list=option_list))
gse <- paste0("GSE",opt$gse)
outFolder <- opt$dir

#===============================================================================
#   Downloading GSE files
#===============================================================================

GEO_dataset <- suppressMessages(getGEO(GEO=gse, GSElimits=NULL, GSEMatrix=TRUE, AnnotGPL=TRUE))

# Iterating into all the GEO datasets
for (i in 1:length(GEO_dataset)) {
	# Exporting expression file as dataframe
	exprs_df <- as.data.frame(exprs(GEO_dataset[[i]]), stringsAsFactors=FALSE)

	# Downloading the annotation file
	annotation <- getGEO(annotation(GEO_dataset[[i]]))

	# Create dictionary probe - gene_name from annotation file
		# Get index for Gene Symbol column
		gene_symbol_col_index <- grep("[Gg]ene.+[Ss]ymbol", colnames(Table(annotation)))

		# Subsetting annotation file and adding probes as rownames
		annotation_df <- as.data.frame(Table(annotation)[,gene_symbol_col_index], stringsAsFactors=FALSE)
		rownames(annotation_df) = Table(annotation)[,1]

	# Binding expression file and annotation df
	expr_df_converted = merge(annotation_df, exprs_df, by="row.names")
	# Formatting final expression file
		expr_df_converted[,1] <- NULL # this removes the probe column
		colnames(expr_df_converted)[1] <- "Gene_name" # renaming the columns containing the gene name
	# Removing duplicated genes, selecting just the most variable
	expr_df_converted <- selectDuplicatedGenes(expr_df_converted)

	# Downloading phenoData
	phenoData <- as.data.frame(pData(GEO_dataset[[i]]), stringsAsFactors=FALSE)
	# Adding Sample Names colums
	phenoData$SampleName <- rownames(phenoData)
	# Renaming characteristics column (most important of the file)
		characteristics_col_index <- grep("[Ss]ource.+[Nn]ame", colnames(phenoData))
		colnames(phenoData)[characteristics_col_index] = "SampleChar"

	#===============================================================================
	#     Exporting files
	#===============================================================================

	# Saving expression matrix
	expr_df_converted.fn = paste0(outFolder,"/eData.",i,".tsv")
	write.table(expr_df_converted, file=expr_df_converted.fn, quote=F, sep="\t", col.names=T, row.names=T)

	# Saving phenoData
	phenoData.fn =  paste0(outFolder,"/pData.",i,".tsv")
	write.table(phenoData, file=phenoData.fn, quote=F, sep="\t", col.names=T, row.names=F)
}
