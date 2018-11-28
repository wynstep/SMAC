###############################################################################
#
#   File name: bcntb.estimate.R
#   Description: Predicts PAM50 breast cancer subtypes using expression matrix
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
suppressMessages(library(plotly))
suppressMessages(library(optparse))
suppressMessages(library(data.table))
source(paste0(GetScriptWD(), "/functions.R")) # file containing functions
source(paste0(GetScriptWD(), "/vars.R")) # file containing vars

#===============================================================================
#    Catching the arguments
#===============================================================================

option_list = list(
  make_option(c("-e", "--exp_file"), action="store", default=NA, type='character',
              help="File containing experimental data"),
  make_option(c("-t", "--target"), action="store", default=NA, type='character',
              help="Clinical data saved in tab-delimited format"),
  make_option(c("-c", "--colouring"), action="store", default=NA, type='character',
              help="Variable from the samples annotation file to be used for samples colouring"),
  make_option(c("-n", "--number"), action="store", default=NA, type='character',
              help="Number of the plot"),
  make_option(c("-d", "--dir"), action="store", default=NA, type='character',
              help="Default directory")
)

opt = parse_args(OptionParser(option_list=option_list))
expFile <- opt$exp_file
annFile <- opt$target
target <- opt$colouring
number <- opt$number
outFolder <- opt$dir
pam50.genes <- fread("src/pam50.genes.txt", sep="\t", header=TRUE, stringsAsFactors= FALSE, data.table=FALSE)

#===============================================================================
#   Pre-processing of the input files
#===============================================================================

  # Read file with expression data
  expData <- read.table(expFile,sep="\t",header=TRUE, stringsAsFactors = FALSE, row.names=NULL)
  rownames(expData) <- make.unique(expData[,1])

  # Read sample annotation file
  annData <- fread(annFile,sep="\t",header=TRUE, data.table=FALSE)

  # Intersecting the samples between expression and target files
  selected_samples <- intersect(as.character(annData[,1]),colnames(expData))

  # subsetting expression file based on selected samples and PAM50 genes (gene name)
  expData.subset <- expData[pam50.genes[,1],colnames(expData) %in% selected_samples]

  # subsetting annData
  annData.subset <- annData[annData[,1] %in% selected_samples,]

  # getting index of target column in the annotation file
  target.index = grep(target, colnames(annData.subset))

  # getting groups names
  groups = unique(annData.subset[,target.index])

  # creating mapping dictionary file_name -- group
  targets.groups = sapply(annData.subset[,1], function(x) {return(annData.subset[annData.subset[,1]==x,target.index])})

#===============================================================================
#     Calculate pam50 subtypes
#===============================================================================

  # creating PAM50 total report dataframe -- empty for now
  pam50.report <- data.frame()

  # iterating all the group present in the Target file
  for (group in groups) {
    # We count the words in the target file, belonging to cancer group
    num_cancer_matches = sum(sapply(CancerToMatch, function(x) grepl(x, group, ignore.case=TRUE)))
    # We repeat the same for the normal terms -- we want to be more stringent, so we'll search just for the exact matches!
    num_normal_matches = sum(sapply(NormalToMatch, function(x) grepl(paste0("\\b",x,"\\b"), group, ignore.case=TRUE)))

    ### Important: if the cancer words are more than normal (Just in that case), we'll perform the analysis
    if (num_cancer_matches > num_normal_matches) {
      # subsetting expression dataframe according to the samples belonging to each group
      expData.group <- as.data.frame(expData.subset[,names(targets.groups[targets.groups==group])])
      rownames(expData.group) = rownames(expData.subset)
      colnames(expData.group) = names(targets.groups[targets.groups==group])

      ## warning, we cannot apply pam50 on a dataframe composed of just one sample
      if (ncol(expData.group) > 1) {
        # Performing PAM50 analysis
        PAM50Preds <- calculatePam50(expData.group, pam50.genes)
        # updating pam50 report
        pam50.report <- rbind(pam50.report, data.frame("File_name"=names(PAM50Preds$subtype),
                                                      "Description"=sapply(names(PAM50Preds$subtype), function(x) targets.groups[[x]]),
                                                      "subtype"=PAM50Preds$subtype))
      }
    }
  }

#===============================================================================
#     Plotting results
#===============================================================================

  ## plots will be produced just if pam50 report is not empty
  if (nrow(pam50.report)>0) {

      ### Barplot ###

        # calculating frequency of different molecular subtypes
        pam50.freq <- as.data.frame(table(pam50.report$subtype))
        # calculating the total number of samples into the report
        total = sum(pam50.freq$Freq)

        # Plotting
        p <- plot_ly(pam50.freq,
                x = pam50.freq$Var1,
                y = pam50.freq$Freq,
                type = 'bar',

              ### defining bar colors ####
              # rgba(255,51,51,0.3) --> red
              # rgba(0,0,204,0.3) --> blue
              # rgba(102,204,0,0.3) --> light green
              # rgba(0,102,0,0.3) --> dark green
              # rgba(204,0,102,0.3) --> dark pink

                marker = list(color = c('rgba(255,51,51,0.3)', 'rgba(0,0,204,0.3)',
                                        'rgba(102,204,0,0.3)', 'rgba(0,102,0,0.3)',
                                        'rgba(204,0,102,0.3)'),
                # defining bar line colors (same as bar colors but different opacity)
                line = list(color = c('rgba(255,51,51,1)', 'rgba(0,0,204,1)',
                                      'rgba(102,204,0,1)', 'rgba(0,102,0,1)',
                                      'rgba(204,0,102,1)'),
                width = 1.5)
              ),
              # adding percentage in the hover label text
              text = paste(round((pam50.freq$Freq/total)*100, digits = 2),'%'),
              textpositionsrc = 'center',
              width = '800px'
        ) %>%
        layout(yaxis = list(title = 'Frequency'), xaxis = list(title = 'Subtypes'), barmode = 'group', title = 'Molecular classification', showlegend = FALSE)

        # Saving as HTML
        widget_fn = paste(outFolder,paste0("molecular_classification.",number,".html"),sep="/")
        saveWidgetFix(p, file=widget_fn, selfcontained = TRUE)

      #===============================================================================
      #     Exporting data
      #===============================================================================

        # saving pam50 results into a text file
        OutFile = paste(outFolder,paste0("target.molecular_classification.",number,".tsv"),sep="/")
        write.table(pam50.report, file=OutFile, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
  }
