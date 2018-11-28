################################################################################
#
#   File name: vars.R
#   Description: This file contains all the important variables for executing BCNTB analyses
#
#   Authors: Stefano Pirro' ( s.pirro@qmul.ac.uk )
#
#   Barts Cancer Institute,
#   Queen Mary, University of London
#   Charterhouse Square, London EC1M 6BQ
#
################################################################################

#===============================================================================
#    Filters for cancers (general) and normal -- collection of Regex for identify cancers
#       the target file
#===============================================================================

CancerToMatch <- c("carcinoma", "cancer", "tumou?r", "metastasis", "sarcoma", "lesion","dysplasia", "atypia", "mastectomy", "immortali.ed")
NormalToMatch <- c("normal", "adjacent", "hyperplasia", "benign", "prophylactic", "mastectomy", "contr.lateral", "control", "wt", "wild.type")
BlacklistToMatch <- c("sample", "rep+", "institute", "university", "centre", "center")

#===============================================================================
#    Receptor genes for assessing the receptor status
#===============================================================================
receptor.genes = c('ESR1','PGR','ERBB2')


#===============================================================================
#    Estimate package
#===============================================================================

  # color bar style
  EstimateCbStyle<-list(
    colorbar = list(title = "Tumour Purity",
                    titleside = "top",
                    titlefont = list(
                      family = "Helvetica Neue",
                      size = 15
                    ),
                    tickfont = list (
                      family = "Helvetica Neue",
                      size = 12,
                      color = 'grey'
                    ),
                    ypad = 10,
                    y = 3,
                    thickness = 15,
                    nticks = 10,
                    xanchor = 'center',
                    xpad = 10),
    size = 5,
    symbol = 'circle',
    cmin=0,
    cmax=1,
    cauto = F
  )

  # margins of the plot
  EstimateMargins <- list(
    l = 10,
    r = 10,
    b = 10,
    t = 10,
    pad = 1
  )
