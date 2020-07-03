################################################################################
#
#   File name: vars.R
#
#   Authors: Stefano Pirro', PhD ( s.pirro@qmul.ac.uk )
#
#   Barts Cancer Institute,
#   Queen Mary, University of London
#   Charterhouse Square, London EC1M 6BQ
#
#		Description: File containing functions for the analysis
#
################################################################################

#===============================================================================
#	Function for cutting too long strings
#===============================================================================
CutString <- function(string, n_words) {
	ul <- unlist(strsplit(string, split="\\s+"))
	if (length(ul) > n_words) {
		paste(ul[1:n_words],collapse=" ")
	} else {
		paste(ul,collapse=" ")
	}
}

#===============================================================================
#	Function for assigning discrete colors to metaData
#===============================================================================
getMetaColours <- function(metaData, typeMap) {
	if (typeMap == "discrete") {
		# Initialising array of pre-defined colours
	  suppressMessages(library(circlize))
	  metaData.colours <- rand_color(100, hue = NULL, luminosity = "random", transparency = 0)
	  # select colours
	  selected.colours = metaData.colours[1:length(unique(metaData))]
	  # creating mapping dictionary colour -- group
	  metaData.colours.map = NULL
	  metaData.colours.map[c(unique(metaData))] <- selected.colours
	  metaData.colours = sapply(metaData, function(x) return(metaData.colours.map[[x]]) )
	  return(metaData.colours)
	}
}

#===============================================================================
#	Function for counting the number of cancer and normal words into an array
#===============================================================================
countCancer <- function(metaData, CancerToMatch, NormalToMatch) {
	# Initialise vector with selected status
	selMD <- c()
	# iterate over the metaData elements and count the cancer/normal terms
	for (mD in metaData) {
		cancerTerms <- sum(sapply(CancerToMatch, function(x) grepl(x, mD, ignore.case=TRUE)))
		normalTerms <- sum(sapply(NormalToMatch, function(x) grepl(x, mD, ignore.case=TRUE)))
		# Append true/false if the row needs to be selected or not
		ifelse (cancerTerms > normalTerms, selMD <- c(selMD, TRUE), selMD <- c(selMD, FALSE))
	}

	# Return the selected rows status
	return(selMD)
}

#===============================================================================
#	Function for determining receptor status
#===============================================================================
calculateRstatus <- function(exprData) {
  # loading required libraries
  suppressMessages(library('mclust'))
  # important step -- filtering genes where the expression value is the same for all the samples, mclust cannot be applied
  exprData <- exprData[!apply(exprData, 1, function(x) all(x==x[1])),]
  # implement MCLUST to a +'ve or -'ve receptor status to each sample, where G=2  i.e. 2-component Gaussian mix
  optimal.model <- apply(exprData, 1, function(x) Mclust(data = x,G = 2, prior = NULL,control = emControl(),initialization = NULL,warn = FALSE))

  # create rstatus report
  rstatus.report <- list()
  for (i in 1:length(rownames(exprData))) {
    gene = rownames(exprData)[i]
    rstatus.report[[gene]] <- optimal.model[[gene]]$classification
  }

	# Transform list of named vectors to dataframe
	rstatus.report.df <- as.data.frame(do.call(rbind, rstatus.report))
	rownames(rstatus.report.df) <- names(rstatus.report)

  # replace 1 values with 0, 2 with 1
  rstatus.report.df[rstatus.report.df == 1] <- 0
  rstatus.report.df[rstatus.report.df == 2] <- 1

  return(rstatus.report.df)
}

#===============================================================================
#	Function for performing molecular classification
#===============================================================================
calculateMclass <- function(exprData, pam50) {
  # loading required libraries
  suppressMessages(library('genefu'))
  # Performing subtyping using PAM50
  PAM50Preds <- molecular.subtyping(sbt.model = "pam50", data=t(exprData), annot=pam50, do.mapping=FALSE)

  return(PAM50Preds)
}

#===============================================================================
#	Function for colouring network gene according to the score
#===============================================================================
ColourNodes <- function (nodes) {
  # loading color palette
  rbPal <- colorRampPalette(c("red", "white", "green"))(length(nodes$color[!is.na(nodes$exprValue)]))

  # colouring selected nodes
  nodes$color[is.na(nodes$exprValue)] = rgb(0.8,0.8,0.8) # light grey
  nodes$color[!is.na(nodes$exprValue)] = rbPal

  # creating label groups
  nodes$group[is.na(nodes$exprValue)] = "Citation NA"
  nodes$group[!is.na(nodes$exprValue)] = ""

  return(nodes)
}

#===============================================================================
# A wrapper to saveWidget which compensates for arguable BUG in
# saveWidget which requires `file` to be in current working
# directory.
#===============================================================================
saveWidgetFix <- function (widget,file,...) {
  # load htmlwidgets library
  suppressMessages(library(htmlwidgets))
	currentWD <- getwd() # store current working directory
	setwd(normalizePath(dirname(file))) # change wd for savewidget
  htmlwidgets::saveWidget(widget,file=basename(file), selfcontained = TRUE,...)
	setwd(currentWD) # Restore working directory
}
