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
#		Description: File containing vars for the analysis
#
################################################################################

#===============================================================================
#   External dependencies
#===============================================================================
dependencies = list(
	"pca" = c("argparse", "data.table", "plotly"),
	"gExpr" = c("argparse", "data.table", "ComplexHeatmap", "circlize"),
	"enrichment" = c("argparse","data.table","clusterProfiler","ReactomePA","rWikiPathways","DOSE","GSEABase","magrittr","msigdbr","ggplot2","enrichplot","org.Hs.eg.db"),
	"rStatus" = c("argparse", "data.table", "plotly"),
	"tPurity" = c("argparse", "data.table", "plotly", "estimate"),
	"mClass" = c("argparse", "data.table", "plotly"),
	"gNet" = c("argparse","data.table","visNetwork","dplyr"),
	"rra" = c("argparse","data.table","RobustRankAggreg","ggplot2"),
	"updateTable" = c("argparse","data.table"),
	"add2mySQL" = c("argparse","data.table","RMariaDB"),
	"dataStats" = c("argparse","data.table","plotly")
)

#===============================================================================
#   Organism libraries
#===============================================================================
organismLibraries <- list(
	"Homo sapiens" = "org.Hs.eg.db", # Human
	"Mus musculus" = "org.Mm.eg.db", # Mouse
	"Rattus norvegicus" = "org.Rn.eg.db", # Rattus
	"Gallus gallus" = "org.Gg.eg.db", # Chicken
	"Caenorhabditis elegans" = "org.Ce.eg.db", # Worm
	"Arabidopsis thaliana" = "org.At.tair.db", # Plant
	"Danio rerio" = "org.Dr.eg.db", # Zebrafish
	"Macaca mulatta" = "org.Mmu.eg.db" # Monkey
)

organismCodes <- list(
	"Homo sapiens" = "hsa", # Human
	"Mus musculus" = "mmu", # Mouse
	"Rattus norvegicus" = "rno", # Rattus
	"Gallus gallus" = "gga", # Chicken
	"Caenorhabditis elegans" = "cel", # Worm
	"Arabidopsis thaliana" = "ath", # Plant
	"Danio rerio" = "dre", # Zebrafish
	"Macaca mulatta" = "mcc" # Monkey
)

#===============================================================================
#   External Commands
#===============================================================================
command = list(
	"img2gif" = "convert -delay 20 -loop 0 $(find %s -name 'update_*.png' -print0 | sort -zV | xargs -r0 echo) %s && rm %s/update_*.png"
)

#===============================================================================
#   Variables
#===============================================================================

allMeshLength <- 29640 # Number of all meSH terms (updated 2020/04/15)
receptorGenes <- c('ESR1','PGR','ERBB2') # Genes for calculating the receptor Status

# cancer/normal keywords to match in the metaData description
CancerToMatch <- c("carcinoma", "cancer", "tumou?r", "metastasis", "sarcoma", "lesion","dysplasia", "atypia", "mastectomy", "immortali.ed")
NormalToMatch <- c("normal", "adjacent", "hyperplasia", "benign", "prophylactic", "mastectomy", "contr.lateral", "control", "wt", "wild.type")

# Margins for plots
plotMargins = list(
	"tPurity" = list( l = 10, r = 10, b = 10, t = 10, pad = 1 )
)

# Style for Plots
EstimateCbStyle<-list(
	colorbar = list(title = "Tumour Purity", titleside = "top",
									titlefont = list( family = "Helvetica Neue", size = 15 ),
									tickfont = list (family = "Helvetica Neue", size = 12, color = 'grey'),
									ypad = 10, y = 3, thickness = 15, nticks = 10, xanchor = 'center', xpad = 10),
	size = 5, symbol = 'circle', cmin=0, cmax=1, cauto = F
)

dbUsername <- "root"
dbPassword <- "hlfy418"
dbName <- "MAP"
