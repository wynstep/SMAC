#! /usr/bin/env python3

##################################################
## SMAC
## Author: Stefano Pirro, PhD (aka wynstep)
## Version: 2.0
## Maintainer: Stefano Pirro, PhD (aka wynstep)
## Email: stefano.pirro@uniroma2.it, s.pirro@mir-nat.com
## Status: Active
##################################################

"""
Description of the script:
	This downloads a list of Papers and relative informations from PubMed.
	It also downloads a series of GEO codes, in order to allow the analysis
	of molecular data in live. This will be performed by the final user, just on
	selected publications.
	PLEASE NOTE: This code is an updated version of SMAC (https://github.com/wynstep/SMAC),
	published on Scientific Reports in 2019 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6642118/)
"""

	## Importing libraries, vars and functions
import os, sys, argparse, datetime, re, gc
import pandas as pd
import numpy as np
from func_timeout import func_timeout, FunctionTimedOut
from vars import *
from functions import *

	## Reading arguments for the analysis
parser = argparse.ArgumentParser()
parser.add_argument('-p','--PMIDlist', help='List of PMIDs to parse', required=True)
parser.add_argument('-a','--analysis', help='analysis to perform (separated by comma) [pca,receptor_status,molecular_classification,tumour_purity,gene_expression,functional_enrichment]', default=defaultAnalyses, required=False)
parser.add_argument('-s','--skipDownload', help='Do you want to skip download of molecular data? <yes,no>', default="no", required=False)
parser.add_argument('-g','--nGenes', help='Number of top differentially-expressed genes to plot in the gene expression analysis', default=100, required=False)
parser.add_argument('-n','--nCores', help='Number of cores for the analysis', default=1, type=int, required=False)
parser.add_argument('-e','--email', help='A valid email address for not being banned by NCBI', required=True)
args = parser.parse_args()

	## Initialising analysis code -- This is the queryTerm without spaces and wrong chars
anCode = os.path.basename(args.PMIDlist).split(".")[0]
currentYear = int(datetime.datetime.now().year) # Logging current year (for MeSH download)
	## Initialising directories
outRootFolder = outRootFolder.format("DATA", anCode)
outGraphDir = outGraphDir.format(outRootFolder)
outGSEDir = outGSEDir.format(outRootFolder)
os.makedirs(outRootFolder, exist_ok = True)
os.makedirs(outGraphDir, exist_ok = True)
os.makedirs(outGSEDir, exist_ok = True)
	## Initialising logFile
logFile = logFile.format(outRootFolder, anCode)

##---------------------------------------##
##	Loading list of PMIDs from argument  ##
##---------------------------------------##
	# We get just unique elements
PMIDs = list(set([element.strip() for element in open(args.PMIDlist,"r").readlines()]))

##--------------------------------------------------------------##
##	Splitting PMIDs into sets for avoiding blocks			    ##
##--------------------------------------------------------------##
PMIDsChunks = list(chunkPMIDs(PMIDs, sizeQuery["chunkPMIDs"]))
	# From this point we do not need the PMIDs variable anymore, so we can delete and save memory
	# This is very useful for low memory machines
gc.collect() # Release garbage memory

##--------------------------------------------------------------##
##	Extracting GEO datasets linked to PMIDs					    ##
##--------------------------------------------------------------##
if args.skipDownload == "no":
		# Initialise file where to store PMID2GSE info
	geoInfoFile = geoInfoFile.format(outRootFolder, anCode)
		# Initialise pandas dataframe that stores data
	PMID2GSE = pd.DataFrame(columns=["PMID","gseCodes","platforms","ftpLinks"])
		# Iterating over each block of PMIDs and retrieve the data (easier for backup purposes)
	for index, pmids in enumerate(PMIDsChunks):
		KeepLog(logFile, "Retrieve GSE codes for batch {0}/{1} \n".format(index, len(PMIDsChunks)))
		# We need to put a timeout to the function, otherwise (if blocked IP or connection, it will stay for forever)
		try:
			tmpPMID2GSE =  func_timeout(timoutThr, RetrieveGeoInfoFromPMID, args=(pmids, logFile))
			SaveIntoDF(tmpPMID2GSE, geoInfoFile, "pandas")
			PMID2GSE = PMID2GSE.append(tmpPMID2GSE, ignore_index=True)
		except FunctionTimedOut:
			KeepLog(logFile, "Sorry, timeout riched for {0}/{1} \n".format(index, len(PMIDsChunks)))
			continue # continue the for loop if function RetrieveGeoInfoFromPMID takes more than 40 second

	# Deleting tmpPMID2GSE, geoInfoFile variable and release garbage memory
	gc.collect()

##--------------------------------------------------------------##
##	Extracting Genes linked to PMIDs						    ##
##--------------------------------------------------------------##
	# Initialise file where to store PMID2Gene info
geneInfoFile = geneInfoFile.format(outRootFolder, anCode)
	# Initialise pandas dataframe that stores data
PMID2Gene = pd.DataFrame(columns=["PMID","geneName","geneScore"])
	# Iterating over each block of PMIDs and retrieve the data (easier for backup purposes)
for index, pmids in enumerate(PMIDsChunks):
	KeepLog(logFile, "Retrieve GENES for batch {0}/{1} \n".format(index, len(PMIDsChunks)))
	# We need to put a timeout to the function, otherwise (if blocked IP or connection, it will stay for forever)
	try:
		tmpPMID2Gene = func_timeout(timoutThr+100, LinkPMID2Gene, args=(pmids, args.email, logFile))
		SaveIntoDF(tmpPMID2Gene, geneInfoFile, "pandas")
		PMID2Gene = PMID2Gene.append(tmpPMID2Gene, ignore_index=True)
	except FunctionTimedOut:
		KeepLog(logFile, "Sorry, timeout riched for {0}/{1} \n".format(index, len(PMIDsChunks)))
		continue # continue the for loop if function RetrieveGeoInfoFromPMID takes more than 40 second

##--------------------------------------------------------------##
##	Create 	PPI (protein-protein interaction network)		    ##
##--------------------------------------------------------------##
if PMID2Gene.shape[0] > 0: # number of rows is more than 0
	PerformAnalysis("gene_network", [geneInfoFile, args.queryTerm, args.nCores, outGraphDir], logFile)
	PerformAnalysis("functional_enrichment_pmid", [geneInfoFile, args.queryTerm, args.nCores, outGraphDir], logFile)

# Deleting tmpPMID2Gene, geneInfoFile variable and release garbage memory
gc.collect()

##--------------------------------------------------------------##
##	Extracting Paper Info									    ##
##--------------------------------------------------------------##

	# Initialising files for storing Papers Info
papersInfoFile = papersInfoFile.format(outRootFolder, anCode)
meshTermsFile = meshTermsFile.format(outRootFolder, anCode)
	# Initialise pandas dataframe and list that stores data
papersInfo = pd.DataFrame(columns=["PMID","author","journal","title","abstract","pubDate","meshHeadings"])
meshInfo = list()

	# Iterating over each block of PMIDs and retrieve the data (easier for backup purposes)
for index, pmids in enumerate(PMIDsChunks):
	KeepLog(logFile, "Retrieve INFO for batch {0}/{1} \n".format(index, len(PMIDsChunks)))
	# We need to put a timeout to the function, otherwise (if blocked IP or connection, it will stay for forever)
	try:
		tmpPapersInfo = func_timeout(timoutThr, RetrievePMIDInfo, args=(pmids, args.email, logFile))
		papersInfo = papersInfo.append(tmpPapersInfo[0], ignore_index=True)
		meshInfo.extend(tmpPapersInfo[1])
		SaveIntoDF(tmpPapersInfo[0], papersInfoFile, "pandas")
		SaveIntoDF(tmpPapersInfo[1], meshTermsFile, "list")
	except FunctionTimedOut:
		KeepLog(logFile, "Sorry, timeout riched for {0}/{1} \n".format(index, len(PMIDsChunks)))
		continue # continue the for loop if function RetrieveGeoInfoFromPMID takes more than 40 second

# Deleting tmpPapersInfo, papersInfoFile variable and release garbage memory
gc.collect()

##--------------------------------------------------------------##
##	Extracting MeSH terms metrics -- If there are NEW mesh term ##
##--------------------------------------------------------------##

if (len(meshInfo) > 0):
	# Initialising files for storing the MeSH statistics
	meshLAboundanceFile = meshLAboundanceFile.format(outRootFolder, anCode)
	meshGAboundanceFile = meshGAboundanceFile.format(outRootFolder, anCode)
	meshLevelFile = meshLevelFile.format(outRootFolder, anCode)
	meshTermsDetailsFile = meshTermsDetailsFile.format(outRootFolder, currentYear, anCode)

	### We first try to download and parse the latest version of MeSH official terms
	### These will be stored into a pandas dataframe
	if (not os.path.isfile(meshTermsDetailsFile)): # Not necessary to parse if the file is present
		officialMeshTerms = DownloadOfficialMesh(currentYear, srcFolder, logFile)
		officialMeshTerms.to_pickle(meshTermsDetailsFile)
	else:
		officialMeshTerms = pd.read_pickle(meshTermsDetailsFile)

	# Load latest version of meshInfo from generated file
	with open(meshTermsFile, 'r') as meshIO:
		meshInfo = meshIO.read().splitlines()
	meshIO.close()

	# Calculation of local aboundance
	lAboundance = RetrieveMeshInfo(meshInfo, "lAboundance", args.email, logFile)
	# Extract JUST official MeSH
	lAboundanceSelected = lAboundance[lAboundance['mesh'].isin(officialMeshTerms['mesh'].tolist())]
	SaveIntoDF(lAboundanceSelected, meshLAboundanceFile, "pandas")

	# Calculating MeSH level for terms present in the paper
	meshLevel = officialMeshTerms[officialMeshTerms['mesh'].isin(lAboundanceSelected['mesh'].tolist())]
	SaveIntoDF(meshLevel, meshLevelFile, "pandas")

	# We'll calculate global aboundance (JUST FOR TERMS WHICH HAVE NOT BEEN CALCULATED)
	gAboundance = pd.DataFrame(columns=['mesh', 'gAboundance'])
	if os.path.isfile(meshGAboundanceFile): # File exists
		meshGAboundanceTerms = pd.read_csv(meshGAboundanceFile, sep="\t")
		# Here we select the terms that have not been previously parsed
		meshGAboundanceTermsSelected = lAboundanceSelected[lAboundanceSelected['mesh'].isin(meshGAboundanceTerms['mesh'].tolist(), invert=True)]
	else: # File does not exists
		meshGAboundanceTermsSelected = lAboundanceSelected

	for j, m in enumerate(meshGAboundanceTermsSelected["mesh"].tolist()):
		KeepLog(logFile, "Retrieve MeSH global aboundance for {0}/{1} -- {2}\n".format(j, len(meshGAboundanceTermsSelected["mesh"].tolist()), m))
		# We need to put a timeout to the function, otherwise (if blocked IP or connection, it will stay for forever)
		try:
			tmpAboundance = func_timeout(timoutThr, RetrieveMeshInfo, args=(m, "gAboundance", args.email, logFile))
			gAboundance.loc[len(gAboundance)] = tmpAboundance
		except FunctionTimedOut:
			KeepLog(logFile, "Sorry, timeout riched for {0}/{1} \n".format(index, len(PMIDsChunks)))
			continue # continue the for loop if function RetrieveGeoInfoFromPMID takes more than 40 second

	SaveIntoDF(gAboundance, meshGAboundanceFile, "pandas")

	##--------------------------------------------------------------##
	##	Perform a Robust Ranking Aggregation on MeSH		 	    ##
	##--------------------------------------------------------------##
	PerformAnalysis("rra", [meshLevelFile, meshGAboundanceFile, meshLAboundanceFile, \
						args.queryTerm, args.nCores, outGraphDir, outRootFolder], logFile)

# Here we clean some variables before running the GSE download and process, since it can be memory expensive
try:
	del meshGAboundanceTerms
	del officialMeshTerms, meshInfo
	del lAboundance, lAboundanceSelected
	del meshLevel
	del gAboundance, meshGAboundanceTermsSelected, tmpAboundance
except:
	KeepLog(logFile, "[WARNING] No need to clear VAR\n")

gc.collect()

##--------------------------------------------------------------##
##	Download GEO data				 							##
##--------------------------------------------------------------##

if args.skipDownload == "no":
		# Iterate over PMID/GSE dataframe and retrieve data
	for index, row in PMID2GSE.iterrows():
		pmid = row.PMID
		gseCodes = [s.strip() for s in row.gseCodes.replace(',', ';').split(';')]
		for gse in gseCodes:
			try:
				gse = "GSE{0}".format(gse)
				gseDir = "{0}/{1}/{2}".format(outGSEDir, pmid, gse)
				KeepLog(logFile, "[MESSAGE] Analysing {0}/{1}\n".format(pmid,gse))
				exprDataFiles, metaDataFile = func_timeout(timoutThr+1000, DownloadGEODataset, args=(gse, gseDir))
			except FunctionTimedOut:
				KeepLog(logFile, "[WARNING] Sorry, timeout riched for {0}/{1}\n".format(pmid, gse))
				subprocess.Popen("rm -rf {0}".format(gseDir), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
				continue # continue the for loop if function RetrieveGeoInfoFromPMID takes more
			except:
				KeepLog(logFile, "[WARNING] Sorry, cannot process {0}/{1}\n".format(pmid, gse))
				subprocess.Popen("rm -rf {0}".format(gseDir), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
				continue

			##--------------------------------------------------------------##
			##	Perform analyses on downloaded datasets				 		##
			##--------------------------------------------------------------##
			if exprDataFiles is not None and metaDataFile is not None: # We perform the analyses JUST if the GEO dataset has been created
				# Iterate over the expression data files generated
				for index, eDataFile in enumerate(exprDataFiles):
					for an in args.analysis.split(","):
						PerformAnalysis(an, [eDataFile, metaDataFile[index], index, args.nCores, args.nGenes, gseDir], logFile)

##--------------------------------------------------------------##
##	Cleaning empty directories				 					##
##--------------------------------------------------------------##
RemoveEmptyDirs("DATA/")
KeepLog(logFile, "SMAC for {0} COMPLETED. Thanks for using it!\n".format(anCode))
