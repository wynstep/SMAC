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

  # subsetting expression file by selected samples and genes (receptor genes for mclust)
  expData.subset <- expData[rownames(expData) %in% receptor.genes,colnames(expData) %in% selected_samples]

  # subsetting annData
  annData.subset <- annData[annData[,1] %in% selected_samples,]

  # getting index of target column in the annotation file
  target.index = grep(target, colnames(annData.subset))

  # getting groups names
  groups = unique(annData.subset[,target.index])

  # creating mapping dictionary file_name -- group
  targets.groups = sapply(annData.subset[,1], function(x) {return(annData.subset[annData.subset[,1]==x,target.index])})

#===============================================================================
#     Analyse receptor status with Mclust
#===============================================================================

  # initilising cont useful for merging expression data
  cont = 1
  # initilising merged expression data
  expData.tumours = data.frame()

  ### Note: the mclust algorithm is applicable just on the tumour samples, not on others
  for (group in groups) {
    # We count the words in the target file, belonging to cancer group
    num_cancer_matches = sum(sapply(CancerToMatch, function(x) grepl(x, group, ignore.case=TRUE)))
    # We repeat the same for the normal terms -- to be more stringent we'll search just for the exact match
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

  # We'll perform the next steps just if there are cancer samples to process and at least 2 receptor genes
  if (ncol(expData.tumours) > 0 & nrow(expData.tumours)>1) {

    # creating Mclust total report dataframe
    mclust.report <- calculateMclust(expData.tumours)
    mclust.report.file <- as.data.frame(t(mclust.report))
    # add sample names and description to final dataframe report
    mclust.report.file$File_name = rownames(mclust.report.file)

    mclust.report.file$Description = sapply(rownames(mclust.report.file), function(x) targets.groups[[x]])
    # write final mclust report
    OutFile = paste(outFolder,paste0("target.receptor_status.",number,".tsv"),sep="/")
    write.table(mclust.report.file, OutFile, sep="\t", row.names = FALSE, col.names=TRUE, quote = FALSE)

    #===============================================================================
    #     Plotting results
    #===============================================================================

      ### Barplot ###

      	# Important note! The triple negative status can be assessed JUST if all the status
      	# for ALL the 3 markers (ER, PR, HER2), is present

        # Preparing dataframe for plotting
          mclust.report.freq <- as.data.frame(apply(mclust.report, 1, function(x) table(x)))

        # Preparing plot filename
          widget_fn = paste(outFolder,paste0("receptor_status.",number,".html"),sep="/")

      	# In all cases, we plot the status for individual markers
        	p.status <- plot_ly(mclust.report.freq,
                          x = colnames(mclust.report.freq),
                          y =  as.numeric(mclust.report.freq[1,]), # first row contains negative values
                          type = 'bar',
                          name = 'Negative',
                          legendgroup="1",
                          # rgba(255,51,51) --> red
                          marker = list(color = 'rgba(255,51,51,0.3)',
                                        line = list(color = 'rgba(255,51,51,1)',
                                                    width = 1.5)),
                          text = paste(round((mclust.report.freq[1,]/sum(mclust.report.freq[,1]))*100, digits = 2),'%'),
                          textpositionsrc = 'center',
                          width = '800px'
                  ) %>%
                  add_trace(mclust.report.freq,
                            x = colnames(mclust.report.freq),
                            y = as.numeric(mclust.report.freq[2,]), # second row contains positive values
                            type = 'bar',
                            name = 'Positive',
                            legendgroup="1",
                            # rgba(0,0,204) --> blue
                            marker = list(color = 'rgba(0,0,204,0.3)',
                                          line = list(color = 'rgba(0,0,204,1)',
                                                      width = 1.5)),
                            text = paste(round((mclust.report.freq[2,]/sum(mclust.report.freq[,1]))*100, digits = 2),'%'),
                            textpositionsrc = 'center',
                            width = '800px'
                  ) %>%
                  layout(yaxis = list(title = 'Count'),
                         barmode = 'stack'
                  )

      	if (nrow(mclust.report) == 3) { # this means we have the status for all the markers...

          # Here we count the number of triple negative samples
          tn.samples = sum(colSums(mclust.report) == 0)
          ntn.samples = sum(colSums(mclust.report) > 0)
          tot = tn.samples + ntn.samples

          p.tneg <- plot_ly(mclust.report.freq,
                                x = "TripleNegative",
                                y = ntn.samples ,
                                type = 'bar',
                                name = 'No',
                                legendgroup="2",
                                 # rgba(102,204,0,0.3) --> light green
                                marker = list(color = 'rgba(102,204,0,0.3)',
                                              line = list(color = 'rgba(102,204,0,1)',
                                                          width = 1.5)),
                                text = paste(round((ntn.samples/tot)*100, digits = 2),'%'),
                                textpositionsrc = 'center',
                                width = '800px'
                        ) %>%
                        add_trace(mclust.report.freq,
                                  x = "TripleNegative",
                                  y = tn.samples,
                                  type = 'bar',
                                  name = 'Yes',
                                  legendgroup="2",
                                  # rgba(204,0,102,0.3) --> dark pink
                                  marker = list(color = 'rgba(204,0,102,0.3)',
                                                line = list(color = 'rgba(204,0,102,1)',
                                                            width = 1.5)),
                                  text = paste(round((tn.samples/tot)*100, digits = 2),'%'),
                                  textpositionsrc = 'center',
                                  width = '800px'
                        ) %>%
                        layout(yaxis = list(title = 'Count'), barmode = 'stack')

          # merging plots together
          mclust.plot <- subplot(p.status, p.tneg, nrows = 1, margin = 0.02, shareY = TRUE)
          mclust.plot <- layout(mclust.plot, title="Receptor Status");

          # Saving as HTML
          saveWidgetFix(mclust.plot, file=widget_fn, selfcontained = TRUE)

      	} else {

          # Saving as HTML
          saveWidgetFix(p.status, file=widget_fn, selfcontained = TRUE)

        }
    }
