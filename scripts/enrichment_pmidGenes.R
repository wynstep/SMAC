################################################################################
#
#   File name: enrichment.R
#
#   Authors: Stefano Pirro', PhD ( s.pirro@qmul.ac.uk )
#
#   Barts Cancer Institute,
#   Queen Mary, University of London
#   Charterhouse Square, London EC1M 6BQ
#
#		Description: this script perform several enrichment analyses, according to different onthologies
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
	source(paste0(GetScriptWD(), "/functions.R")) # file containing functions
	source(paste0(GetScriptWD(), "/vars.R")) # file containing vars
	# Taking advantage of the pacman module for installing unmet dependencies
	if (!require("pacman")) install.packages("pacman")
	pacman::p_load(dependencies[["enrichment"]], character.only = T)
	org <- "Homo sapiens"

#===============================================================================
#	Load Arguments for the analysis
#===============================================================================
	parser <- ArgumentParser()
	parser$add_argument("-i", "--geneFile", action="store", required=TRUE, dest="geneFile", help="Input gene list to analyse")
	parser$add_argument("-l", "--queryLabel", action="store", required=TRUE, dest="queryLabel", help="Label to assign to output file")
	parser$add_argument("-n", "--nCores", action="store", required=TRUE, dest="nCores", type="integer", default=1, help="Number of available cores for the analysis")
	parser$add_argument("-o", "--outputDir", action="store", required=TRUE, dest="outputDir", help="Output directory")
	args <- parser$parse_args()

#===============================================================================
#	Load input files
#===============================================================================
	# load gene list
	geneData <- fread(args$geneFile, sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE, fill=TRUE, data.table=FALSE, nThread=args$nCores)

#===============================================================================
#	Convert Gene names to Entrez ID
#===============================================================================
	geneConverted <- bitr(geneData$geneName, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
	geneList <- geneConverted$ENTREZID

#===============================================================================
#	WikiPathways enrichment
#===============================================================================
	# Logging operation
	cat("[MESSAGE] WikiPathways enrichment...\n")
	# create folder for WikiPathway analysis
	wpDir <- paste(args$outputDir,"WikiPathway",sep="/")
	dir.create(wpDir, showWarnings = FALSE)
	tryCatch({
		# We first download the pathway data associated to the selected organism
		pathFile <- downloadPathwayArchive(organism=org, format="gmt", destpath=wpDir)
		# Loading pathFile and parsing with clusterProfiler
		wpgmtfile <- paste(wpDir, pathFile, sep="/")
		wp2gene <- read.gmt(wpgmtfile)
		# take advantage of dplyr and tidyr for extracting Pathways enriched for the selected genes
		wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
		wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
		wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME
		ewp <- enricher(geneList, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
		# Converting entrez ids to gene names
		ewp <- setReadable(ewp, organismLibraries[["Homo sapiens"]], keyType = "ENTREZID")
		# Save WikiPathway enrichment into a file
		ewpFile <- paste(wpDir,paste("WikiPathway","enrichment","tsv",sep="."), sep="/")
		write.table(ewp, file=ewpFile, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
		# Generate Plots
			# DotPlot
		dotPlot <- dotplot(ewp, showCategory=50) + ggtitle("DotPlot for WikiPathways (top 50)")
		dotPlotFN <- paste(wpDir, paste("WikiPathway","dotplot","png", sep="."), sep="/")
		ggsave(dotPlot, file=dotPlotFN, dpi=600, width=50, height=50, unit="cm")
			# HeatMap
		heatmap <- heatplot(ewp) + ggtitle("Heatmap for WikiPathways")
		heatmapFN <- paste(wpDir, paste("WikiPathway","heatmap","png", sep="."), sep="/")
		ggsave(heatmap, file=heatmapFN, dpi=600, width=80, height=40, unit="cm", scale=0.8)
			# UpSet plot
		upset <- upsetplot(ewp) + ggtitle("UpSet plot for WikiPathways")
		upsetFN <- paste(wpDir, paste("WikiPathway","upset","png", sep="."), sep="/")
		ggsave(upset, file=upsetFN, dpi=600, width=80, height=40, unit="cm", scale=0.8)
			# Delete temporary file
		file.remove(wpgmtfile)
	}, error = function(e) {
		cat("[ERROR] Cannot perform WikiPathways enrichment\n")
		system(paste("rm -rf",wpDir,sep=" "))
	})

#===============================================================================
#	MSigDb enrichment
#===============================================================================
	# Logging operation
	cat("[MESSAGE] MSigDb enrichment...\n")
	# create folder for MSigDb analysis
	msDir <- paste(args$outputDir,"MSigDb",sep="/")
	dir.create(msDir, showWarnings = FALSE)
	tryCatch({
		# Download gene set associated to the organism and subset columns
		m_df <- msigdbr(species = org, category = "C5") %>% dplyr::select(gs_name, entrez_gene)
		# Perform enrichment
		em <- enricher(geneList, TERM2GENE=m_df)
		setReadable(em, org.Hs.eg.db, keyType = "ENTREZID")
		# Save MSigDB enrichment into a file
		emFile <- paste(msDir,paste("MSigDB","enrichment","tsv",sep="."), sep="/")
		write.table(em, file=emFile, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
		# Generate Plots
			# DotPlot
		dotPlot <- dotplot(em, showCategory=50) + ggtitle("DotPlot for MSigDB (top 50)")
		dotPlotFN <- paste(msDir, paste("MSigDB","dotplot","png", sep="."), sep="/")
		ggsave(dotPlot, file=dotPlotFN, dpi=600, width=50, height=50, unit="cm")
			# HeatMap
		heatmap <- heatplot(em) + ggtitle("Heatmap for MSigDB")
		heatmapFN <- paste(msDir, paste("MSigDB","heatmap","png", sep="."), sep="/")
		ggsave(heatmap, file=heatmapFN, dpi=600, width=80, height=40, unit="cm", scale=0.8)
			# UpSet plot
		upset <- upsetplot(em) + ggtitle("UpSet plot for MSigDB")
		upsetFN <- paste(msDir, paste("MSigDB","upset","png", sep="."), sep="/")
		ggsave(upset, file=upsetFN, dpi=600, width=80, height=40, unit="cm", scale=0.8)
	}, error = function(e) {
		cat("[ERROR] Cannot perform MSigDb enrichment\n")
		# print(e)
		system(paste("rm -rf",msDir,sep=" "))
	})

#===============================================================================
#	Disease onthology enrichment
#===============================================================================
	# Logging operation
	cat("[MESSAGE] DiseaseOnthology enrichment...\n")
	# create folder for DiseaseOnthology analysis
	doDir <- paste(args$outputDir,"DiseaseOnthology",sep="/")
	dir.create(doDir, showWarnings = FALSE)
	tryCatch({
		# Perform the enrichment
		edo <- enrichDO(gene = geneList, ont = "DO", pvalueCutoff = 0.05,
										pAdjustMethod = "BH", minGSSize = 5, maxGSSize = 500,
										qvalueCutoff = 0.05, readable = TRUE)
		# Save DO enrichment into a file
		edoFile <- paste(doDir,paste("DiseaseOnthology","enrichment","tsv",sep="."), sep="/")
		write.table(edo, file=edoFile, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
		# Generate Plots
			# DotPlot
		dotPlot <- dotplot(edo, showCategory=50) + ggtitle("DotPlot for DiseaseOnthology (top 50)")
		dotPlotFN <- paste(doDir, paste("DiseaseOnthology","dotplot","png", sep="."), sep="/")
		ggsave(dotPlot, file=dotPlotFN, dpi=600, width=50, height=50, unit="cm")
			# HeatMap
		heatmap <- heatplot(edo) + ggtitle("Heatmap for DiseaseOnthology")
		heatmapFN <- paste(doDir, paste("DiseaseOnthology","heatmap","png", sep="."), sep="/")
		ggsave(heatmap, file=heatmapFN, dpi=600, width=80, height=40, unit="cm", scale=0.8)
			# UpSet plot
		upset <- upsetplot(edo) + ggtitle("UpSet plot for DiseaseOnthology")
		upsetFN <- paste(doDir, paste("DiseaseOnthology","upset","png", sep="."), sep="/")
		ggsave(upset, file=upsetFN, dpi=600, width=80, height=40, unit="cm", scale=0.8)
	}, error = function(e) {
		cat("[ERROR] Cannot perform DiseaseOnthology enrichment\n")
		# print(e)
		system(paste("rm -rf",doDir,sep=" "))
	})

#===============================================================================
#	Network of Cancer Gene enrichment
#===============================================================================
	# Logging operation
	cat("[MESSAGE] NetworkCancerGene enrichment...\n")
	# create folder for DiseaseOnthology analysis
	ncgDir <- paste(args$outputDir,"NetworkCancerGene",sep="/")
	dir.create(ncgDir, showWarnings = FALSE)
	tryCatch({
		# Perform the enrichment
		ncg <- enrichNCG(geneList, readable = TRUE)
		# Save NCG enrichment into a file
		ncgFile <- paste(ncgDir,paste("NetworkCancerGene","enrichment","tsv",sep="."), sep="/")
		write.table(ncg, file=ncgFile, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
		# Generate Plots
			# DotPlot
		dotPlot <- dotplot(ncg, showCategory=50) + ggtitle("DotPlot for NetworkCancerGene (top 50)")
		dotPlotFN <- paste(ncgDir, paste("NetworkCancerGene","dotplot","png", sep="."), sep="/")
		ggsave(dotPlot, file=dotPlotFN, dpi=600, width=50, height=50, unit="cm")
			# HeatMap
		heatmap <- heatplot(ncg) + ggtitle("Heatmap for NetworkCancerGene")
		heatmapFN <- paste(ncgDir, paste("NetworkCancerGene","heatmap","png", sep="."), sep="/")
		ggsave(heatmap, file=heatmapFN, dpi=600, width=80, height=40, unit="cm", scale=0.8)
			# UpSet plot
		upset <- upsetplot(ncg) + ggtitle("UpSet plot for NetworkCancerGene")
		upsetFN <- paste(ncgDir, paste("NetworkCancerGene","upset","png", sep="."), sep="/")
		ggsave(upset, file=upsetFN, dpi=600, width=80, height=40, unit="cm", scale=0.8)
	}, error = function(e) {
		cat("[ERROR] Cannot perform NetworkCancerGene enrichment\n")
		# print(e)
		system(paste("rm -rf",ncgDir,sep=" "))
	})

#===============================================================================
#	DisGeNET enrichment
#===============================================================================
	# Logging operation
	cat("[MESSAGE] DisGeNET enrichment...\n")
	# create folder for DiseaseOnthology analysis
	dgnDir <- paste(args$outputDir,"DisGeNET",sep="/")
	dir.create(dgnDir, showWarnings = FALSE)
	tryCatch({
		# Perform the enrichment
		dgn <- enrichDGN(geneList, readable=TRUE)
		# Save NCG enrichment into a file
		dgnFile <- paste(dgnDir,paste("DisGeNET","enrichment","tsv",sep="."), sep="/")
		write.table(dgn, file=dgnFile, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
		# Generate Plots
			# DotPlot
		dotPlot <- dotplot(dgn, showCategory=50) + ggtitle("DotPlot for DisGeNET (top 50)")
		dotPlotFN <- paste(dgnDir, paste("DisGeNET","dotplot","png", sep="."), sep="/")
		ggsave(dotPlot, file=dotPlotFN, dpi=600, width=50, height=50, unit="cm")
			# HeatMap
		heatmap <- heatplot(dgn) + ggtitle("Heatmap for DisGeNET")
		heatmapFN <- paste(dgnDir, paste("DisGeNET","heatmap","png", sep="."), sep="/")
		ggsave(heatmap, file=heatmapFN, dpi=600, width=80, height=40, unit="cm", scale=0.8)
			# UpSet plot
		upset <- upsetplot(dgn) + ggtitle("UpSet plot for DisGeNET")
		upsetFN <- paste(dgnDir, paste("DisGeNET","upset","png", sep="."), sep="/")
		ggsave(upset, file=upsetFN, dpi=600, width=80, height=40, unit="cm", scale=0.8)
	}, error = function(e) {
		cat("[ERROR] Cannot perform DisGeNET enrichment\n")
		# print(e)
		system(paste("rm -rf",dgnDir,sep=" "))
	})

#===============================================================================
#	Gene onthology enrichment
#===============================================================================
	# create folder for DiseaseOnthology analysis
	goDir <- paste(args$outputDir,"GeneOnthology",sep="/")
	dir.create(goDir, showWarnings = FALSE)
	tryCatch({
		# Logging operation
		cat("[MESSAGE] GeneOnthologyCC enrichment...\n")
		# Perform the enrichment
		ego.cc <- enrichGO(gene = geneList, OrgDb = organismLibraries[[org]],
											ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.01,
											qvalueCutoff = 0.05, readable = TRUE) # Cellular Compartment (CC)
		# Save EGO CC enrichment into a file
		egoccFile <- paste(goDir,paste("GeneOnthologyCC","enrichment","tsv",sep="."), sep="/")
		write.table(ego.cc, file=egoccFile, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
		# Generate Plots
			# DotPlot
		dotPlot <- dotplot(ego.cc, showCategory=50) + ggtitle("DotPlot for GeneOnthologyCC (top 50)")
		dotPlotFN <- paste(goDir, paste("GeneOnthologyCC","dotplot","png", sep="."), sep="/")
		ggsave(dotPlot, file=dotPlotFN, dpi=600, width=50, height=50, unit="cm")
			# HeatMap
		heatmap <- heatplot(ego.cc) + ggtitle("Heatmap for GeneOnthologyCC")
		heatmapFN <- paste(goDir, paste("GeneOnthologyCC","heatmap","png", sep="."), sep="/")
		ggsave(heatmap, file=heatmapFN, dpi=600, width=80, height=40, unit="cm", scale=0.8)
			# UpSet plot
		upset <- upsetplot(ego.cc) + ggtitle("UpSet plot for GeneOnthologyCC")
		upsetFN <- paste(goDir, paste("GeneOnthologyCC","upset","png", sep="."), sep="/")
		ggsave(upset, file=upsetFN, dpi=600, width=80, height=40, unit="cm", scale=0.8)

		# Logging operation
		cat("[MESSAGE] GeneOnthologyBP enrichment...\n")
		ego.bp <- enrichGO(gene = geneList, OrgDb = organismLibraries[[org]],
											ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01,
											qvalueCutoff = 0.05, readable = TRUE) # Biological Processes (BP)
		# Save EGO BP enrichment into a file
		egobpFile <- paste(goDir,paste("GeneOnthologyBP","enrichment","tsv",sep="."), sep="/")
		write.table(ego.bp, file=egobpFile, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
		# Generate Plots
			# DotPlot
		dotPlot <- dotplot(ego.bp, showCategory=50) + ggtitle("DotPlot for GeneOnthologyBP (top 50)")
		dotPlotFN <- paste(goDir, paste("GeneOnthologyBP","dotplot","png", sep="."), sep="/")
		ggsave(dotPlot, file=dotPlotFN, dpi=600, width=50, height=50, unit="cm")
			# HeatMap
		heatmap <- heatplot(ego.bp) + ggtitle("Heatmap for GeneOnthologyBP")
		heatmapFN <- paste(goDir, paste("GeneOnthologyBP","heatmap","png", sep="."), sep="/")
		ggsave(heatmap, file=heatmapFN, dpi=600, width=80, height=40, unit="cm", scale=0.8)
			# UpSet plot
		upset <- upsetplot(ego.bp) + ggtitle("UpSet plot for GeneOnthologyBP")
		upsetFN <- paste(goDir, paste("GeneOnthologyBP","upset","png", sep="."), sep="/")
		ggsave(upset, file=upsetFN, dpi=600, width=80, height=40, unit="cm", scale=0.8)

		# Logging operation
		cat("[MESSAGE] GeneOnthologyMF enrichment...\n")
		ego.mf <- enrichGO(gene = geneList, OrgDb = organismLibraries[[org]],
											ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.01,
											qvalueCutoff = 0.05, readable = TRUE) # Molecular Function (MF)
		# Save EGO MF enrichment into a file
		egomfFile <- paste(goDir,paste("GeneOnthologyMF","enrichment","tsv",sep="."), sep="/")
		write.table(ego.mf, file=egomfFile, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
		# Generate Plots
			# DotPlot
		dotPlot <- dotplot(ego.mf, showCategory=50) + ggtitle("DotPlot for GeneOnthologyMF (top 50)")
		dotPlotFN <- paste(goDir, paste("GeneOnthologyMF","dotplot","png", sep="."), sep="/")
		ggsave(dotPlot, file=dotPlotFN, dpi=600, width=50, height=50, unit="cm")
			# HeatMap
		heatmap <- heatplot(ego.mf) + ggtitle("Heatmap for GeneOnthologyMF")
		heatmapFN <- paste(goDir, paste("GeneOnthologyMF","heatmap","png", sep="."), sep="/")
		ggsave(heatmap, file=heatmapFN, dpi=600, width=80, height=40, unit="cm", scale=0.8)
			# UpSet plot
		upset <- upsetplot(ego.mf) + ggtitle("UpSet plot for GeneOnthologyMF")
		upsetFN <- paste(goDir, paste("GeneOnthologyMF","upset","png", sep="."), sep="/")
		ggsave(upset, file=upsetFN, dpi=600, width=80, height=40, unit="cm", scale=0.8)
	}, error = function(e) {
		cat("[ERROR] Cannot perform GeneOnthology enrichment\n")
		# print(e)
		system(paste("rm -rf",goDir,sep=" "))
	})

#===============================================================================
#	KEGG pathway enrichment
#===============================================================================
	# create folder for DiseaseOnthology analysis
	keggDir <- paste(args$outputDir,"KEGG",sep="/")
	dir.create(keggDir, showWarnings = FALSE)
	tryCatch({
		# Logging operation
		cat("[MESSAGE] KEGG enrichment...\n")
		# Perform the enrichment
		kk <- enrichKEGG(gene = geneList, organism = organismCodes[[org]], pvalueCutoff = 0.05)
		# Save KEGG enrichment into a file
		kkFile <- paste(keggDir,paste("KEGG","enrichment","tsv",sep="."), sep="/")
		write.table(kk, file=kkFile, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
		# Generate Plots
			# DotPlot
		dotPlot <- dotplot(kk, showCategory=50) + ggtitle("DotPlot for KEGG (top 50)")
		dotPlotFN <- paste(keggDir, paste("KEGG","dotplot","png", sep="."), sep="/")
		ggsave(dotPlot, file=dotPlotFN, dpi=600, width=50, height=50, unit="cm")
			# HeatMap
		heatmap <- heatplot(kk) + ggtitle("Heatmap for KEGG")
		heatmapFN <- paste(keggDir, paste("KEGG","heatmap","png", sep="."), sep="/")
		ggsave(heatmap, file=heatmapFN, dpi=600, width=80, height=40, unit="cm", scale=0.8)
			# UpSet plot
		upset <- upsetplot(kk) + ggtitle("UpSet plot for KEGG")
		upsetFN <- paste(keggDir, paste("KEGG","upset","png", sep="."), sep="/")
		ggsave(upset, file=upsetFN, dpi=600, width=80, height=40, unit="cm", scale=0.8)
	}, error = function(e) {
		cat("[ERROR] Cannot perform KEGG enrichment\n")
		# print(e)
		system(paste("rm -rf",keggDir,sep=" "))
	})

#===============================================================================
#	Reactome pathway enrichment
#===============================================================================
	# create folder for DiseaseOnthology analysis
	reactDir <- paste(args$outputDir,"REACTOME",sep="/")
	dir.create(reactDir, showWarnings = FALSE)
	tryCatch({
		# Logging operation
		cat("[MESSAGE] REACTOME enrichment...\n")
		# Reactome enrichment
		rpa <- enrichPathway(gene = geneList, pvalueCutoff=0.05)
		# Save REACTOME enrichment into a file
		rpaFile <- paste(reactDir,paste("REACTOME","enrichment","tsv",sep="."), sep="/")
		write.table(rpa, file=rpaFile, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
		# Generate Plots
			# DotPlot
		dotPlot <- dotplot(rpa, showCategory=50) + ggtitle("DotPlot for REACTOME (top 50)")
		dotPlotFN <- paste(reactDir, paste("REACTOME","dotplot","png", sep="."), sep="/")
		ggsave(dotPlot, file=dotPlotFN, dpi=600, width=50, height=50, unit="cm")
			# HeatMap
		heatmap <- heatplot(rpa) + ggtitle("Heatmap for REACTOME")
		heatmapFN <- paste(reactDir, paste("REACTOME","heatmap","png", sep="."), sep="/")
		ggsave(heatmap, file=heatmapFN, dpi=600, width=80, height=40, unit="cm", scale=0.8)
			# UpSet plot
		upset <- upsetplot(rpa) + ggtitle("UpSet plot for REACTOME")
		upsetFN <- paste(reactDir, paste("REACTOME","upset","png", sep="."), sep="/")
		ggsave(upset, file=upsetFN, dpi=600, width=80, height=40, unit="cm", scale=0.8)
	}, error = function(e) {
		cat("[ERROR] Cannot perform REACTOME enrichment\n")
		# print(e)
		system(paste("rm -rf",reactDir,sep=" "))
	})
