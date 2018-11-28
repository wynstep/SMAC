#title           :CreateNetwork.R
#description     :Create interaction network for the genes correlated to the literature
#author          :s.pirro
#date            :20180123
#version         :0.1
#notes           :
#==============================================================================

# loading packages
suppressMessages(library(optparse))
suppressMessages(library(visNetwork))
suppressMessages(library(data.table))

# setup arguments
option_list <- list(
  make_option(c("-l", "--list"),
              help="lists of genes"),
  make_option(c("-n", "--network"),
              help="network file (mentha format)"),
  make_option(c("-o", "--output"),
              help="output directory")
)
opt <- parse_args(OptionParser(option_list=option_list))

geneFile = opt$list
outFolder = opt$output
netFile = opt$network

#===============================================================================
#   Custom functions
#===============================================================================

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
  rbPal <- colorRampPalette(c("red", "white", "green"))(length(nodes$color[!is.na(nodes$expr_value)]))

  # colouring selected nodes
  nodes$color[is.na(nodes$expr_value)] = rgb(0.8,0.8,0.8) # light grey
  nodes$color[!is.na(nodes$expr_value)] = rbPal

  # creating label groups
  nodes$group[is.na(nodes$expr_value)] = "Citation NA"
  nodes$group[!is.na(nodes$expr_value)] = ""

  return(nodes)
}

#===============================================================================
#   Pre-processing of the input files
#===============================================================================

  #=============================================================================
  #     Genes File
  #=============================================================================

    # Read gene file
    geneData <- fread(geneFile,sep="\t",header=FALSE, stringsAsFactors = FALSE, data.table=FALSE, fill=TRUE)
    colnames(geneData) = c("Gene", "Score")
    # making gene names unique
    geneData$Gene = make.unique(geneData$Gene)

  #=============================================================================
  #     Network File
  #=============================================================================

    # Read Network file
    netData <- fread(netFile,sep=";",header=TRUE, stringsAsFactors = FALSE, data.table=FALSE)
    # Filtering Network file according to selected genes and thresholds
    netData.filtered <- netData[(netData[,2] %in% geneData[,1] | netData[,5] %in% geneData[,1]) & (netData[,7] >= 0.7 & netData[,7] <= 1), ]

    # Collecting total list of genes selected
    all.genes <- data.frame("Gene" = unique(c(netData.filtered[,2], netData.filtered[,5])), stringsAsFactors=FALSE)

#===============================================================================
#   Building network
#===============================================================================

  # mapping gene file with all the genes in the network
  mapped.genes = merge(all.genes, geneData, by="Gene", all.x=T)
  # check if the merge succeed, otherwise treat genData as mapped
  if (nrow(mapped.genes) == 0) {
    mapped.genes = geneData
  }

  #================================================
  # Creating Network data
  #================================================

  # creating genes (nodes) dataframe -- for visNetwork
  nodes <- data.frame(
    id=mapped.genes$Gene,
    expr_value = mapped.genes$Score,
    label=mapped.genes$Gene,
    color="black", # just a sample color
    group="gene",
    shape="circle",
    title = paste0("<a href='http://www.genecards.org/cgi-bin/carddisp.pl?gene=",mapped.genes$Gene,"' target=new><b>", mapped.genes$Gene ,"</b></a><br>
                   <b> Enrichment score (PubMed): </b>",round(mapped.genes$Score,digits=3),""),
    shadow = FALSE,
    stringsAsFactors = FALSE)

  # Coloring genes (nodes) according to the pubmed value
  nodes.coloured <- ColourNodes(nodes, genes)

  # creating edges -- for visNetwork
  edges <- data.frame(from = netData.filtered[,2], 
                      to = netData.filtered[,5], 
                      length = netData.filtered[,7], 
                      title = paste0("<b> Source: </b>",netData.filtered[,2]," <b> Target: </b>", netData.filtered[,5] ,"<br>
                                      <b> Interaction score: </b>",round(netData.filtered[,7],digits=3),"")
                      )
  #================================================
  # Saving Network
  #================================================
  WidgetFile = paste0(outFolder,"/gene_network.html")

  visNetwork(nodes.coloured, edges, height = "600px", width = "100%") %>%
  visEdges(smooth = TRUE, physics = TRUE) %>%
  visPhysics(stabilization = FALSE, solver = "barnesHut") %>%
  # save newtwork into file
  saveWidgetFix(file = WidgetFile, background = "white", selfcontained = TRUE)
