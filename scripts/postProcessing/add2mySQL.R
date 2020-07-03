################################################################################
#
#   File name: add2mySQL.R
#
#   Authors: Stefano Pirro', PhD ( s.pirro@mir-nat.com)
#
#
#		Description: Starting from the papersInfo table, this script
#   updates the corrensponding mySQL table
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
	pacman::p_load(dependencies[["add2mySQL"]], character.only = T)
  options(useFancyQuotes = FALSE) # Removing fancy quotes from quoting

#===============================================================================
#	Load Arguments for the analysis
#===============================================================================
	parser <- ArgumentParser()
	parser$add_argument("-a", "--papersFile", action="store", required=TRUE, help="TSV containing papers details")
  parser$add_argument("-t", "--table", action="store", required=TRUE, help="Name of the table to update")
  args <- parser$parse_args()

#===============================================================================
#	Load input data
#===============================================================================
  papersData <- fread(args$papersFile, sep="\t", header=TRUE, stringsAsFactors=FALSE, fill=TRUE, data.table=FALSE)

#===============================================================================
#	Establish connection with database
#===============================================================================
  papersDBcon <- dbConnect(RMariaDB::MariaDB(), username=dbUsername, password=dbPassword, dbname=dbName, host="localhost")

#===============================================================================
#	Running query for adding data
#===============================================================================

  for (i in 1:nrow(papersData)) {
    # Initialise variables
    PMID <- sQuote(gsub("'", "''", papersData[i,"PMID"]))
    author <- sQuote(gsub("'", "''", papersData[i,"author"]))
    journal <- sQuote(gsub("'", "''", papersData[i,"journal"]))
    title <- sQuote(gsub("'", "''", papersData[i,"title"]))
    abstract <- sQuote(gsub("'", "''", papersData[i,"abstract"]))
    pubDate <- sQuote(gsub("'", "''", papersData[i,"pubDate"]))
    meshHeadings <- sQuote(gsub("'", "''", papersData[i,"meshHeadings"]))
    geneName <- sQuote(gsub("'", "''", papersData[i,"geneName"]))
    availableAnalyses <- sQuote(gsub("'", "''", papersData[i,"availableAnalyses"]))
    # Initialise query
    query <- paste0("INSERT IGNORE INTO ",args$table," (PMID,author,journal,title,abstract,pubDate,meshHeadings,geneName,availableAnalyses) VALUES (",
                    paste(PMID,author,journal,title,abstract,pubDate,meshHeadings,geneName,availableAnalyses,sep=",")
                    ,");"
                  )
    cat("Adding row ",i,"/",nrow(papersData),"to ",args$table,"\n")
    dbInsert <- dbSendQuery(papersDBcon, query)
    # Clear the result.
    dbClearResult(dbInsert)
  }

	# Trasform pubDate column to Date
	dateQuery <- paste0("UPDATE ",args$table," SET pubDate = str_to_date( pubDate, '%d/%m/%Y' )")
	dbInsert <- dbSendQuery(papersDBcon, dateQuery)
	# Clear the result.
	dbClearResult(dbInsert)

  # Disconnect to clean up the connection to the database.
  dbDisconnect(papersDBcon)
