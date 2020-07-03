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
import os, sys, argparse, datetime, itertools, re, gc, glob
import pandas as pd
import numpy as np
from func_timeout import func_timeout, FunctionTimedOut
from vars import *
from functions import *

	## Reading arguments for the analysis
parser = argparse.ArgumentParser()
parser.add_argument('-t','--geneList', help='List of genes for starting a PubMed search query', required=True)
parser.add_argument('-l','--PMIDsLimit', help='Max number of PMIDs to retrieve (default = all)', type=int, default=None, required=False)
parser.add_argument('-a','--analysis', help='analysis to perform (separated by comma) [pca,receptor_status,molecular_classification,tumour_purity,gene_expression,functional_enrichment]', default=defaultAnalyses, required=False)
parser.add_argument('-s','--skipDownload', help='Do you want to skip download of molecular data? <yes,no>', default="no", required=False)
parser.add_argument('-g','--nGenes', help='Number of top differentially-expressed genes to plot in the gene expression analysis', default=100, required=False)
parser.add_argument('-n','--nCores', help='Number of cores for the analysis', default=1, type=int, required=False)
parser.add_argument('-e','--email', help='A valid email address for not being banned by NCBI', required=True)
args = parser.parse_args()

	## Initialising analysis code -- This is the queryTerm without spaces and wrong chars
anCode = os.path.basename(args.geneList).split(".")[0]
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
##	Converting genes to nucleotide codes ##
##---------------------------------------##
	# Initialise pandas dataframe that stores data
GENE2NUCL = pd.DataFrame(columns=["GeneID","GeneName","NucleotideID"])
	# Initialise file for storing gene conversion
nuclInfoFile = nuclInfoFile.format(outRootFolder, anCode)
	# Try to recover analysis if already done
if (os.path.isfile(nuclInfoFile)):
	# Load previous gene set
	with open(nuclInfoFile, "r") as nuclInfoData:
		for line in nuclInfoData:
			line = line.rstrip().split("\t")
			pNuclInfo.append("{0}\t{1}".format(line[0],line[1]))
	# Load current gene set
	cNuclInfo = []
	with open(args.geneList, "r") as geneData:
		for line in geneData:
			line = line.rstrip().split("\t")
			cNuclInfo.append("{0}\t{1}".format(line[0],line[1]))
	# Calculate the intersection between current and previous Nucl Info
	newGenes = list(set(cNuclInfo) - set(pNuclInfo))
	KeepLog(logFile, "[MESSAGE] Recovery mode ON --> Old Genes: {0} New Genes: {1}\n".format(len(pNuclInfo), len(newGenes)))
else:
	newGenes = []
	with open(args.geneList, "r") as geneData:
		for line in geneData:
			line = line.rstrip().split("\t")
			newGenes.append("{0}\t{1}".format(line[0],line[1]))

	# Iterating over Gene info and retrieve data
for index, geneInfo in enumerate(newGenes):
	geneID = geneInfo.rstrip().split("\t")[0]
	geneName = geneInfo.rstrip().split("\t")[1]
	KeepLog(logFile, "Retrieve Nucleotide code for batch {0}/{1}\n".format(index, len(newGenes)))
	# We need to put a timeout to the function, otherwise (if blocked IP or connection, it will stay for forever)
	try:
		tmpGENE2NUCL =  func_timeout(timoutThr, Gene2Nucleotide, args=(geneID, geneName, args.email))
		if tmpGENE2NUCL is not None:
			SaveIntoDF(tmpGENE2NUCL, nuclInfoFile, "pandas")
		else:
			continue
	except FunctionTimedOut:
		KeepLog(logFile, "Sorry, timeout riched for {0}/{1} \n".format(index, len(newGenes)))
		continue # continue the for loop if function RetrieveGeoInfoFromPMID takes more than 40 second

##---------------------------------------##
##	Connecting Nucleotide to PMIDs 		 ##
##---------------------------------------##
	# Initialise pandas dataframe that stores data
NUCL2PMID = pd.DataFrame(columns=["GeneID","GeneName","NucleotideID","PMID"])
	# Initialise file for storing gene conversion
nuclPmidFile = nuclPmidFile.format(outRootFolder, anCode)
	# Search for presence of nuclInfoFile
if (os.path.isfile(nuclInfoFile)):
	nuclInfoData = pd.read_csv(nuclInfoFile, sep="\t")
	# Check if nuclPmidFile has been already generated (backup purposes)
	if (os.path.isfile(nuclPmidFile)):
		KeepLog(logFile, "[MESSAGE] Recovery mode ON\n")
		nuclPmidData = pd.read_csv(nuclPmidFile, sep="\t")
		pNuclID = nuclPmidData["NucleotideID"].tolist()
	# Iterating over Nucleotide codes and retrieve PMIDs correlated with the gene
	for index, row in nuclInfoData.iterrows():
		geneID = row["GeneID"]
		geneName = row["GeneName"]
		nucleotideID = row["NucleotideID"]
		# Check if current nucleotideID has been already analysed
		if nucleotideID not in pNuclID:
			# We need to put a timeout to the function, otherwise (if blocked IP or connection, it will stay for forever)
			try:
				KeepLog(logFile, "Retrieve PMIDs for {0}/{1} \n".format(index, nuclInfoData.shape[0]))
				tmpNucl2PMID =  func_timeout(timoutThr, Nucleotide2PMID, args=(geneID, geneName, nucleotideID, args.email))
				SaveIntoDF(tmpNucl2PMID, nuclPmidFile, "pandas")
				NUCL2PMID = NUCL2PMID.append(tmpNucl2PMID, ignore_index=True)
			except FunctionTimedOut:
				KeepLog(logFile, "Sorry, timeout riched for {0}/{1} \n".format(index, nuclInfoData.shape[0]))
				continue # continue the for loop if function RetrieveGeoInfoFromPMID takes more than 40 second
			except:
				continue
		else:
			KeepLog(logFile, "{0}/{1} already analysed\n".format(index, nuclInfoData.shape[0]))

##---------------------------------------##
##	Aggregate PMIDs and Genes	 		 ##
##---------------------------------------##
	# Check presence of nuclPmidFile
if (os.path.isfile(nuclPmidFile)):
	# Load into pandas dataframe
	nuclPmidData = pd.read_csv(nuclPmidFile, sep="\t")
	# Initialise file for storing aggregated PMID conversion
	nuclPmidFile_aggregated = nuclPmidFile_aggregated.format(outRootFolder, anCode)
	# Group by PMID and aggregate genes
	nuclPmidData['NucleotideID'] = nuclPmidData['NucleotideID'].astype(str)
	nuclPmidData_aggregated = nuclPmidData.groupby(["PMID"], as_index = False).agg({'GeneID': ','.join, 'GeneName': ','.join, 'NucleotideID': ','.join})
	SaveIntoDF(nuclPmidData_aggregated, nuclPmidFile_aggregated, "pandas")
	# Take not of all PMIDs to analyse
	PMIDs = nuclPmidData_aggregated["PMID"].tolist()

##--------------------------------------------------------------##
##	Splitting PMIDs into sets for avoiding blocks			    ##
##--------------------------------------------------------------##
PMIDsChunks = list(chunkPMIDs(PMIDs, sizeQuery["chunkPMIDs"]))
	# From this point we do not need the PMIDs variable anymore, so we can delete and save memory
	# This is very useful for low memory machines
gc.collect() # Release garbage memory

##---------------------------------------##
##	Download informations from PMIDs	 ##
##---------------------------------------##
	# Iterating over each block of PMIDs and retrieve the data (easier for backup purposes)
for index, pmids in enumerate(PMIDsChunks):
	# We need to put a timeout to the function, otherwise (if blocked IP or connection, it will stay for forever)
	try:
		KeepLog(logFile, "Retrieve INFO for batch {0}/{1} \n".format(index, len(PMIDsChunks)))
		tmpPapersInfo = func_timeout(timoutThr, RetrievePMIDInfo, args=(pmids, args.email, logFile))
		# Save papers info and MeSH term into file
		SaveIntoDF(tmpPapersInfo[0], papersInfoFile, "pandas")
		SaveIntoDF(tmpPapersInfo[1], meshTermsFile, "list")
		# Update dataframes
		papersInfo = papersInfo.append(tmpPapersInfo[0], ignore_index=True)
		meshInfo.extend(tmpPapersInfo[1])
	except FunctionTimedOut:
		KeepLog(logFile, "Sorry, timeout riched for {0}/{1} \n".format(index, len(PMIDsChunks)))
		continue # continue the for loop if function RetrieveGeoInfoFromPMID takes more than 40 second
	except Exception as e:
		#print(e)
		continue

# Deleting tmpPapersInfo variable and release garbage memory
gc.collect()

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
