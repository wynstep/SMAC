################################################################################
#
#   File name: functions.R
#   Description: This file contains all the important functions for executing BCNTB analyses
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

source(paste0(GetScriptWD(), "/vars.R")) # file containing vars

# This function classifies breast cancer samples, according to the pam50 signature
calculatePam50 <- function(expData, pam50.genes) {

  # loading required libraries
  suppressMessages(library('genefu'))

  # Performing subtyping using PAM50
  PAM50Preds <- molecular.subtyping(sbt.model = "pam50", data=t(expData), annot=pam50.genes, do.mapping=FALSE)

  return(PAM50Preds)
}

# This function applies the mclust algorithm to the selected expression dataframe
calculateMclust <- function(expData) {
  # loading required libraries
  suppressMessages(library('mclust'))

  # important step -- filtering genes where the expression value is the same for all the samples, mclust cannot be applied
  expData <- expData[!apply(expData, 1, function(x) all(x==x[1])),]
  # implement MCLUST to a +'ve or -'ve receptor status to each sample, where G=2  i.e. 2-component Gaussian mix
  optimal.model <- apply(expData, 1, function(x) Mclust(data = x,G = 2, prior = NULL,control = emControl(),initialization = NULL,warn = FALSE))

  # create mclust report
  mclust.report <- list()
  for (i in 1:length(rownames(expData))) {
    gene = rownames(expData)[i]
    mclust.report[[i]] <- as.data.frame(optimal.model[[gene]]$classification)
  }

  # transform list into dataframe, then transpose
  mclust.report <- t(data.frame(mclust.report))
  rownames(mclust.report) = rownames(expData)

  # replace 1 values with 0, 2 with 1
  mclust.report[mclust.report == 1] <- 0
  mclust.report[mclust.report == 2] <- 1

  return(mclust.report)
}

# Assign colours to analysed groups
getTargetsColours <- function(targets) {
  # Initialising array of pre-defined colours
  suppressMessages(library(randomcoloR))
  targets.colours <- distinctColorPalette(k=100, altCol = FALSE, runTsne = FALSE)

  # select colours
  selected.colours = targets.colours[1:length(unique(targets))]

  # creating mapping dictionary colour -- group
  targets.colours.map = NULL
  targets.colours.map[c(unique(targets))] <- selected.colours

  # mapping targets elements to colours
  targets.colours = sapply(targets, function(x) {return(targets.colours.map[[x]])})

  # returning values
  return(targets.colours)
}

# A wrapper to saveWidget which compensates for arguable BUG in
# saveWidget which requires `file` to be in current working
# directory.
saveWidgetFix <- function (widget,file,...) {
  # load htmlwidgets library
  suppressMessages(library(htmlwidgets))
  wd<-getwd()
  on.exit(setwd(wd))
  outDir<-dirname(file)
  file<-basename(file)
  setwd(outDir);
  htmlwidgets::saveWidget(widget,file=file,...)
}

# Function for colouring nodes according to the zscore
ColourNodes <- function (nodes, genes) {
  # loading color palette
  rbPalpos <- colorRampPalette(c(rgb(1,0.8,0.4),rgb(1,0,0)))(length(nodes$color[!is.na(nodes$expr_value) & nodes$expr_value > 0]))
  rbPalneg <- colorRampPalette(c(rgb(0,0,1),rgb(0,0.7,1)))(length(nodes$color[!is.na(nodes$expr_value) & nodes$expr_value <= 0]))

  # colouring selected nodes
  nodes$color[is.na(nodes$expr_value)] = rgb(0.8,0.8,0.8) # light grey
  nodes$color[!is.na(nodes$expr_value) & nodes$expr_value > 0] = rbPalpos
  nodes$color[!is.na(nodes$expr_value) & nodes$expr_value <= 0] = rbPalneg
  nodes$color[nodes$label %in% genes] = "yellow" #rgb(0.8,1,0) # yellow

  # creating label groups
  nodes$group[is.na(nodes$expr_value)] = "Expression NA"
  nodes$group[!is.na(nodes$expr_value) & nodes$expr_value > 0] = "Up-regulated"
  nodes$group[!is.na(nodes$expr_value) & nodes$expr_value <= 0] = "Down-regulated"
  nodes$group[nodes$label %in% genes] = "Selection"

  return(nodes)
}

# This function computes the "information score" for each element in the pData matrix
computeInfoScore <- function(info) {
  # parse info string as character
  info <- as.character(info)
  # splitting the info string by space
  info.vector <- strsplit(info, "\\s+")
  # calculating the number of words composing the string
  info.length = length(info.vector[[1]])
  # Counting the number of Cancer-releated words in the string
  num_cancer_matches = sum(sapply(CancerToMatch, function(x) grepl(x, info, ignore.case=TRUE)))
  # We repeat the same for the normal terms -- to be more stringent we'll search for the exact match
  num_normal_matches = sum(sapply(NormalToMatch, function(x) grepl(paste0("\\b",x,"\\b"), info, ignore.case=TRUE)))
  # We repeat the same for the blacklisted terms
  num_blacklist_matches = sum(sapply(BlacklistToMatch, function(x) grepl(x, info, ignore.case=TRUE)))

  ## Calculating the scores
  if (num_blacklist_matches > 0) {
    score = -1
  } else {
    score =  (num_cancer_matches+num_normal_matches)/info.length
  }

  # checking presence of NaN and replace with 0
  score = ifelse(is.na(score), 0, score)
  return(score)
}
