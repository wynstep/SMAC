################################################################################
#
#   File name: bcntb.estimate.R
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
suppressMessages(library(estimate))
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

#===============================================================================
#   Pre-processing of the input files
#===============================================================================

  # Read file with expression data
  expData <- read.table(expFile,sep="\t",header=TRUE, stringsAsFactors = FALSE, row.names=NULL)
  rownames(expData) = make.unique(expData[,1])

  # Read sample annotation file
  annData <- fread(annFile,sep="\t",header=TRUE, data.table=FALSE)

  # Intersecting the samples between expression and target files
  selected_samples <- intersect(as.character(annData[,1]),colnames(expData))

  # subsetting expression file by selected samples
  expData.subset <- expData[,colnames(expData) %in% selected_samples]

  # subsetting annData
  annData.subset <- annData[annData[,1] %in% selected_samples,]

  # getting index of target column in the annotation file
  target.index = grep(target, colnames(annData.subset))

  # getting groups names
  groups = unique(annData.subset[,target.index])

  # creating mapping dictionary file_name -- group
  targets.groups = sapply(annData.subset[,1], function(x) {return(annData.subset[annData.subset[,1]==x,target.index])})


#===============================================================================
#     Estimate tumour purity
#===============================================================================

  # initilising cont useful for merging expression data
  cont = 1
  # initilising merged expression data
  expData.tumours = data.frame()

  ### Note: the mclust algorithm is applicable just on the tumour samples, not on others
  for (group in groups) {
    # We count the words in the target file, belonging to cancer group
    num_cancer_matches = sum(sapply(CancerToMatch, function(x) grepl(x, group, ignore.case=TRUE)))
    # We repeat the same for the normal terms -- here we want to be more stringent, so we'll search for the exact match
    num_normal_matches = sum(sapply(NormalToMatch, function(x) grepl(paste0("\\b",x,"\\b"), group, ignore.case=TRUE)))
    ### Important: if the cancer words are more than normal (Just in that case), we'll perform the analysis
    if (num_cancer_matches > num_normal_matches) {
      # filtering the expression data
      expData.tmp = as.data.frame(expData.subset[,names(targets.groups[targets.groups==group])])
      rownames(expData.tmp) = rownames(expData.subset)
      colnames(expData.tmp) = names(targets.groups[targets.groups==group])

      if (cont == 1) {
        expData.tumours = expData.tmp
      } else {
        expData.tumours = cbind(expData.tumours, expData.tmp)
      }
      # increment cont
      cont = cont +1
    }
  }

  # We'll perform the next steps just if there are cancer samples to process
  if (nrow(expData.tumours) > 0) {

      # the "scale" function performs z-score scaling by columns, for this reason we had to transpose (t) the matrix,
      # perform the scaling, then transpose again
      expData.tumours.scaled <- as.data.frame(t(scale(t(expData.tumours))))

      # writing expression data filtered on a tmp file (needed by estimate package)
      expData.tmp.fn = paste(outFolder,'expData.txt',sep="/")
      write.table(expData.tumours.scaled,file=expData.tmp.fn,row.names=T,col.names = T,sep="\t",quote=F)

      # filtering common genes using the built-in Estimate function
      estFile = paste(outFolder,"CommonGenesEstimate.gct",sep="/")
      suppressMessages(filterCommonGenes(expData.tmp.fn, estFile, id = "GeneSymbol"))

      # estimate purity score
      scoreFile = paste(outFolder,"EstimateScore.gct",sep="/")
      suppressMessages(estimateScore(estFile, output.ds=scoreFile))

    #===============================================================================
    #     Plotting results
    #===============================================================================

      ### Scatterplot ###

        # reading the file generated by Estimate and import into a dataframe
        estData=read.table(file=scoreFile,sep="\t",skip=2,header=T,row.names=1)
        estData.df = data.frame('File_name'=colnames(estData[,-1]),
                                'Description'=sapply(colnames(estData[,-1]), function(x) targets.groups[[x]]),
                                'StromalScore'=round(as.numeric(estData[1,-1]), digits = 4),
                                'ImmuneScore'=round(as.numeric(estData[2,-1]), digits = 4),
                                'ESTIMATEScore'=round(as.numeric(estData[3,-1]), digits = 4),
                                'TumourPurity'=round(as.numeric(estData[4,-1]), digits = 4)
                                )

        # Plotting
        p <- plot_ly( estData.df,
                      x= ~StromalScore,
                      y= ~ImmuneScore,
                      z= ~ESTIMATEScore,
                      color = ~TumourPurity,
                      hoverinfo= 'text',
                      text = paste('Sample Name:', estData.df$File_name,
                         '<br><br>StromalScore:', estData.df$StromalScore,
                         '<br>ImmuneScore:', estData.df$ImmuneScore,
                         '<br>ESTIMATEScore:', estData.df$ESTIMATEScore,
                         '<br><br><b>Tumour Purity</b>:', estData.df$TumourPurity),
                      marker = EstimateCbStyle,
                      width = '800px' ) %>%
              add_markers() %>%
              layout(scene = list(xaxis = list(title = 'StromalScore', orientation = "v"),
                                  yaxis = list(title = 'ImmuneScore'),
                                  zaxis = list(title = 'EstimateScore')),
                    margin = EstimateMargins,
                    title = "")

        # Saving as HTML
        widget_fn = paste(outFolder,paste0("tumour_purity.",number,".html"),sep="/")
        saveWidgetFix(p, file=widget_fn, selfcontained = TRUE)

    #===============================================================================
    #     Exporting data
    #===============================================================================

      # saving Estimate results into a text file
      OutFile = paste(outFolder,paste0("target.tumour_purity.",number,".tsv"),sep="/")
      write.table(estData.df, file=OutFile, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

    #===============================================================================
    #     Cleaning residues
    #===============================================================================

    # removing gct files generated by the estimate package
    gct_files <- dir(path=outFolder, pattern="gct")
    file.remove(paste0(outFolder, "/", gct_files))
    # removing tmp expression file
    file.remove(expData.tmp.fn)
  }
