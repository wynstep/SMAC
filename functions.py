#! /usr/bin/env python3

##################################################
## SMAC
## File containing the functions for the analysis

## Author: Stefano Pirro, PhD (aka wynstep)
## Version: 2.0
## Maintainer: Stefano Pirro, PhD (aka wynstep)
## Email: stefano.pirro@uniroma2.it, s.pirro@mir-nat.com
## Status: Active
##################################################

## Importing libraries, vars and functions
import os, sys, GEOparse, re, shutil, subprocess, json
import pandas as pd
import numpy as np
import xml.etree.ElementTree as et # Library for parsing XML files
from tenacity import *
from dateutil import parser # library to convert any type of text into datetime format
from Bio import Entrez, Medline
from vars import *

##--------------------------------------------------------------##
##	Skipping in case of ERROR								    ##
##--------------------------------------------------------------##
def return_last_value(retry_state):
	pass

##--------------------------------------------------------------##
##	Print on terminal 										    ##
##--------------------------------------------------------------##
def terminalPrint(message):
	sys.stdout.write(message)
	sys.stdout.flush()

##--------------------------------------------------------------##
##	Keep Log of all the operations							    ##
##--------------------------------------------------------------##
def KeepLog(logFile, message):
	# print log message on the terminal first
	terminalPrint(message)
	# open a log stream
	logIO = open(logFile, "a")
	# write the log message
	logIO.write(message)
	# close stream
	logIO.close()

##--------------------------------------------------------------##
##	Save extracted datasets into dataframes (TSV files)		    ##
##--------------------------------------------------------------##
def SaveIntoDF(data, outFile, type):
	if type == "pandas":
		if not os.path.isfile(outFile):
			data.to_csv(outFile, sep='\t', index=False, encoding='utf-8')
		else: # else it exists so append without writing the header
			data.to_csv(outFile, sep='\t', mode='a', index=False, header=False, encoding='utf-8')

	elif type == "list":
		outIO = open(outFile, "a")
		for d in data:
			outIO.write("{0}\n".format(d))
		outIO.close()

##--------------------------------------------------------------##
##	Check if folder exists and is empty		    				##
##--------------------------------------------------------------##
def checkDir(directory):
	# Check if directory exists
	if os.path.exists(directory) and not os.path.isfile(directory):
		# Checking if the directory is empty or not
		if not os.listdir(directory):
			nFiles = 0
		else:
			# Count number of Files into the directory
			nFiles = len(os.listdir(directory))
	else:
		nFiles = 0

	return nFiles

##--------------------------------------------------------------##
##	Hard cleaning garbage memory			    				##
##--------------------------------------------------------------##
def cleanMemory():
	cleanCmd = "echo \"echo 3 > /proc/sys/vm/drop_caches\""
	subprocess.check_call(cleanCmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

##--------------------------------------------------------------##
##	Remove empty directories								    ##
##--------------------------------------------------------------##
@retry(stop=stop_after_attempt(7), wait=wait_random(min=5, max=7), retry_error_callback=return_last_value)
def RemoveEmptyDirs(rootDir):
	# Iterating over elements in the rootDir
	for root, dirs, files in os.walk(rootDir, topdown=False):
		for dir in dirs:
			if len(os.listdir( os.path.join(root, dir) )) == 0: #check whether the directory is empty
				os.rmdir( os.path.join(root, dir) )

##--------------------------------------------------------------##
##	Chunk list of PMIDs into chunks of N size				 	##
##--------------------------------------------------------------##
def chunkPMIDs(pmids, n):
	# For item i in a range that is a length of l,
	for i in range(0, len(pmids), n):
		# Create an index range for l of n items:
		yield pmids[i:i+n]

##--------------------------------------------------------------##
##	Link gene names to nucleotide codes				 			##
##--------------------------------------------------------------##
@retry(stop=stop_after_attempt(2), wait=wait_random(min=5, max=7), retry_error_callback=return_last_value)
def Gene2Nucleotide(geneID, geneName, email):
	# initialise email for not being blocked by NCBI
	Entrez.email = email
	# Take advantage of bioPython e-utils for retrieving the list of nucleotide codes
	# associated with gene names
	totalNC = Entrez.esearch(db='nucleotide', term=geneID, retmax=2)
	NCData = Entrez.read(totalNC)
	# Initialise tmpGENE2NUCL
	tmpGENE2NUCL = pd.DataFrame([[geneID, geneName, NCData["IdList"][0]]], columns=["GeneID","GeneName","NucleotideID"])
	# Return dataframe
	return tmpGENE2NUCL

##--------------------------------------------------------------##
##	Link nucleotide IDs to PMIDs		 						##
##--------------------------------------------------------------##
@retry(stop=stop_after_attempt(2), wait=wait_random(min=5, max=7), retry_error_callback=return_last_value)
def Nucleotide2PMID(geneID, geneName, nucleotideID, email):
	# initialise email for not being blocked by NCBI
	Entrez.email = email
	# Initialise empty pandas dataframe
	tmpNUCL2PMID = pd.DataFrame(columns=["GeneID","GeneName","NucleotideID","PMID"])
	# Link nucleotide to PMID
	pmidQuery = Entrez.elink(dbfrom='nucleotide', db='pubmed', id=nucleotideID, cmd='neighbor_score')
	pmidInfo = Entrez.read(pmidQuery)
	# Iterate over the results and extract the data
	index = 0
	for pi in pmidInfo:
		if (len(pi["LinkSetDb"]) > 0): # this means we found some results
			for pmid in pi["LinkSetDb"][0]["Link"]:
				tmpNUCL2PMID.loc[index] = [geneID, geneName, nucleotideID, pmid["Id"]]
				index = index + 1 # Update location

	return tmpNUCL2PMID



##--------------------------------------------------------------##
##	Retrieve Papers from PubMed								    ##
##--------------------------------------------------------------##
@retry(stop=stop_after_attempt(7), wait=wait_random(min=5, max=7), retry_error_callback=return_last_value)
def RetrievePMIDs(term, limit, email, logFile):
	# The email is very important for not being blocked by PubMed in the retrieval process
	Entrez.email = email
	# Clearing the container of PMIDs (declared in the vars.py)
	pubmedIDs.clear()
	# Take advantage of bioPython e-utils for retrieving the list of PMIDs
	# Papers are sorted according to the "relevance", measure defined by NCBI authors
	# https://pubmed.ncbi.nlm.nih.gov/30153250-best-match-new-relevance-search-for-pubmed/
	totalPubs = Entrez.esearch(db='pubmed', term=term, sort="relevance")
	# Reading data retrieved from PubMed
	PubsData = Entrez.read(totalPubs)
	# Logging Number of retrieved publications
	logMessage = "Total number of publications: {0} \n".format(PubsData["Count"])
	KeepLog(logFile, logMessage)

	## Since the e-utilities server allows to extract a maximum of 10,000 PMIDs per query,
	## We need to split the number of retrieved items in chunks

	# check if limit is a number or not
	if limit == None:
		for i in range(0,int(PubsData["Count"]),sizeQuery["RetrievePMIDs"]):
			KeepLog(logFile, "Retrieve PMIDs for {0}/{1} \n".format(i, int(PubsData["Count"])))
			totalPubs = Entrez.esearch(db='pubmed', retmax=sizeQuery["RetrievePMIDs"], retstart=i, term=term)
			# reading results
			PubsData = Entrez.read(totalPubs)
			# Save list of PMIDs into the pubmedIDs array
			pubmedIDs.extend(PubsData['IdList'])
	else:
		size = sizeQuery["RetrievePMIDs"] if limit > sizeQuery["RetrievePMIDs"] else limit
		for i in range(0,limit,size):
			KeepLog(logFile, "Retrieve PMIDs for {0}/{1} \n".format(i, limit))
			totalPubs = Entrez.esearch(db='pubmed', retmax=size, retstart=i, term=term)
			# reading results
			PubsData = Entrez.read(totalPubs)
			# Save list of PMIDs into the pubmedIDs array
			pubmedIDs.extend(PubsData['IdList'])

	return pubmedIDs

##--------------------------------------------------------------##
##	Retrieve GEO datasets IDs from each PMID code				##
##--------------------------------------------------------------##
@retry(stop=stop_after_attempt(2), wait=wait_random(min=5, max=7), retry_error_callback=return_last_value)
def RetrieveGeoInfoFromPMID(PMIDs, email):
	# The email is very important for not being blocked by PubMed in the retrieval process
	Entrez.email = email

	# Initialise tmpPMID2GSE
	tmpPMID2GSE = pd.DataFrame(columns=["PMID","gseCodes","platforms","ftpLinks"])
	# Entrez query
	geoQuery = Entrez.elink(dbfrom='pubmed', db='gds', id=PMIDs)
	geoInfo = Entrez.read(geoQuery)

	# Iterate over the results and extract the data
	start = 0
	for gi in geoInfo:
		gdsCodes.clear() # Initialise container with GDS codes
		pmid = gi["IdList"][0]
		if (len(gi["LinkSetDb"]) > 0): # this means we found some results
			for gds in gi["LinkSetDb"][0]["Link"]:
				gdsCodes.append(gds["Id"])

			# retrieve GSE codes and other info from GDS
			geoQuery = Entrez.esummary(db='gds', id=",".join(list(set(gdsCodes))))
			geoInfo = Entrez.read(geoQuery)
			geoQuery.close()

			# initialising arrays
			gseCodes, platforms, ftpLinks = ([] for i in range(3))

			# iterating retrieved output and extract infos
			for i in range(0, len(geoInfo)):
				gseCodes.append(geoInfo[i].get("GSE", "NA"))
				platforms.append(geoInfo[i].get("GPL", "NA"))
				ftpLinks.append(geoInfo[i].get("FTPLink", "NA"))

			# uniquing arrays and returning
			tmpPMID2GSE.loc[start] = [pmid, \
					','.join(str(x) for x in list(set(gseCodes))),\
					','.join(str(x) for x in list(set(platforms))),\
					','.join(str(x) for x in list(set(ftpLinks)))
				]
			# incresing index for updating at different locations
			start = start + 1

	return tmpPMID2GSE

##--------------------------------------------------------------##
##	Retrieve gene info (single set of PMIDs)				    ##
##--------------------------------------------------------------##
@retry(stop=stop_after_attempt(2), wait=wait_random(min=2, max=3)) #, retry_error_callback=return_last_value)
def LinkPMID2Gene(PMIDs, email, logFile):
	# Initialise query to Entrez
	Entrez.email = email
	geneQuery = Entrez.elink(dbfrom='pubmed', db='gene', id=PMIDs, cmd='neighbor_score')
	geneInfo = Entrez.read(geneQuery)
	geneScores = {} # Initialise container with gene scores

	# Iterate over the results and extract the data
	for gi in geneInfo:
		# keep track of PMID
		pmid = gi["IdList"][0]

		if (len(gi["LinkSetDb"]) > 0): # this means we found some results
			for gene in gi["LinkSetDb"][0]["Link"]:
				if int(gene["Score"]) > 0:
					geneScores[gene["Id"]] = (int(gene["Score"]), pmid)

	# Retrieve detailed informations from gene list
	geneDetQuery = Entrez.esummary(db='gene', id=",".join(list(geneScores.keys())), version="2.0", retmode="json")
	geneDetInfo = json.load(geneDetQuery)
	tmpPMID2Gene = pd.DataFrame(columns=["PMID","geneName","geneScore"]) # initialise pandas dataframe

	# Parse retrieved information and add to pandas dataframe
	start = 0
	for gene, scorePMID in geneScores.items():
		try:
			tmpPMID2Gene.loc[start] = [scorePMID[1], geneDetInfo["result"][gene]["name"], scorePMID[0]]
			start = start + 1
		except:
			next

	return tmpPMID2Gene

##--------------------------------------------------------------##
##	Retrieve papers info (single set of PMIDs)				    ##
##--------------------------------------------------------------##
#@retry(stop=stop_after_attempt(2), wait=wait_random(min=5, max=7), retry_error_callback=return_last_value)
def RetrievePMIDInfo(PubmedIDs, email, logFile):

	# We'll take advantage of ePost function for speeding up the process
	Entrez.email = email
	searchHandle = Entrez.epost(db='pubmed', id=",".join(str(pmid) for pmid in PubmedIDs))
	searchResults = Entrez.read(searchHandle)
	# Getting web environment and query_key (for efetch after)
	webenv, queryKey = searchResults["WebEnv"], searchResults["QueryKey"]

	# Clearing the container of Mesh Terms (declared in the vars.py)
	meshHeadingsContainer.clear()
	# Initialising pandas dataframe with papersInfo
	tmpPapersInfo = pd.DataFrame(columns=["PMID","author","journal","title","abstract","pubDate","meshHeadings"])
	# Initialising cont for updating pandas dataframe
	cont = 0

	# Data will be retrieved in chucks of PMIDs of fixed length
	PMIDsInfo = Entrez.efetch(db='pubmed', rettype='medline', retmode='text', retmax=len(PubmedIDs), webenv=webenv, query_key=queryKey)
	records = Medline.parse(PMIDsInfo)

	## iterating into records set
	for record in records:
		# Getting vars and assign defualt if not exist
		pmid = record.get('PMID', '')
		au = ",".join(record.get('AU', '')) #author
		journal = record.get('JT', '') # journal title
		title = record.get('TI', '') #title
		abstract = record.get('AB', '') # abstract
		meshHeadings = [mh.replace('*', '') for mh in record.get('MH', '')] # mesh terms

		# in case of strage dates, we add 1st January 2020 as default
		try:
			pubDate = parser.parse(record.get('DP', '')).strftime("%Y-%m-%d")
		except ValueError:
			pubDate = parser.parse("2020-01-01").strftime("%Y-%m-%d")

		# updating pandas dataframe
		tmpPapersInfo.loc[cont]=[pmid,au,journal,title,abstract,pubDate,",".join(meshHeadings)]
		# updating mesh heading list
		meshHeadingsContainer.extend(meshHeadings)
		# Update cont
		cont = cont + 1

	# Return variables
	return tmpPapersInfo, meshHeadingsContainer

##--------------------------------------------------------------##
##	Retrieve meSH terms info			    					##
##--------------------------------------------------------------##
@retry(stop=stop_after_attempt(1), wait=wait_random(min=5, max=7), retry_error_callback=return_last_value)
def RetrieveMeshInfo(meshList, type, email, logFile):

	uniqueMesh, counts = np.unique(meshList, return_counts=True)
	lAboundanceDict = dict(zip(uniqueMesh, counts))

	if type == "lAboundance":
		##--------------------------------------------------------------##
		##	Get local Aboundance				    					##
		##--------------------------------------------------------------##
		KeepLog(logFile, "Retrieve MeSH local aboundance \n")
		# Create pandas dataframe from dict
		aboundance = pd.DataFrame(lAboundanceDict.items(), columns=['mesh', 'lAboundance'])

	elif type == "gAboundance":

		##--------------------------------------------------------------##
		##	Get global Aboundance				    					##
		##--------------------------------------------------------------##

		## Since this is a longer and more risky process, we call the function
		# for each single mesh term.

			# Initialising email for Entrez queries
		Entrez.email = email
			# Querying Entrez database
		query = Entrez.esearch(db='pubmed', term="{0} [MeSH Terms]".format(meshList), rettype='Count')
		result = Entrez.read(query)
			# store results in a pandas row
		aboundance = [meshList, result['Count']]

	return aboundance

##--------------------------------------------------------------##
##	Download GEO data											##
##--------------------------------------------------------------##
@retry(stop=stop_after_attempt(2), wait=wait_random(min=1, max=2), retry_error_callback=return_last_value)
def DownloadGEODataset(gse, gseDir):
	# Query for retriecing GEO data
	gseData = GEOparse.get_GEO(geo=gse, destdir=gseDir, how='full', annotate_gpl=True, include_data=True, silent=True)
	# Initialise data containers
	exprData = pd.DataFrame(columns=["ID_REF"])
	metaData = pd.DataFrame()
	exprDataMapped = []
	exprDataFiles = []
	metaDataFiles = []

	# Iterating over all GSM samples in GSE and collect data JUST if RNA-seq
	for gsmName, gsm in gseData.gsms.items():
		if gsm.metadata['type'][0]=='RNA' and len(gsm.table)>0: # We are looking at RNA-seq data with stored data
			##----------------------##
			##	Expression			##
			##----------------------##
				# extract expression data into a pandas dataframe
			tmpExpr = pd.DataFrame({ "ID_REF":list(gsm.table['ID_REF']), gsmName:list(gsm.table['VALUE']) })
				# Appending to exprData (by column)
			exprData = tmpExpr if exprData.shape[0] == 0 else exprData.merge(tmpExpr, how='outer', on='ID_REF')

			##----------------------##
			##	Metadata			##
			##----------------------##
				# extract metadata into a pandas dataframe
			tmpMetadata = pd.DataFrame(gsm.metadata.items(), columns=['MetaData', gsmName])
				# manipulate metaData for being used by subsequent analyses
			tmpMetadata = tmpMetadata.set_index('MetaData') # make metadata column as row name
			tmpMetadata = tmpMetadata.transpose()
				# Appending to metaData (by column)
			metaData = tmpMetadata if metaData.shape[0] == 0 else metaData.append(tmpMetadata)

		# We'll map the platform data JUST if we have data collected in the exprData
	if exprData.shape[0] > 0:
		##----------------------##
		##	platform			##
		##----------------------##
			# Here we aim at mapping the ID_REF to the gene name
		for gplName, gpl in gseData.gpls.items():
				# Get name of column containing gene symbol and ID information
			idCol = [col for col in gpl.table.columns if re.search("^ID", col, re.IGNORECASE)]
			gsCol = [col for col in gpl.table.columns if re.search("^gene(.+|)(symbol|name|id)", col, re.IGNORECASE)]
				# Check if gsCol is present
			if len(gsCol) > 0 and len(idCol) > 0:
				platformData = pd.DataFrame({ "ID_REF": list(gpl.table[idCol[0]]), "geneName": list(gpl.table[gsCol[0]]), })
				##----------------------------------##
				##	Map genes to gene platforms		##
				##----------------------------------##
				exprDataMapped.append(pd.merge(platformData, exprData, how="inner", on="ID_REF"))
			else:
				exprDataMapped.append(exprData)

	##-----------------------------------------------##
	##	Save datasets into files -- or delete folder ##
	##-----------------------------------------------##

		# Iterate over mapped expression data
	if len(exprDataMapped) > 0:
		for index, eData in enumerate(exprDataMapped):
			if eData.shape[0] > 0 and metaData.shape[0] > 0:
				exprDataFile = "{0}/{1}.exprs.{2}.tsv".format(gseDir, gse, index)
				metaDataFile = "{0}/{1}.meta.{2}.tsv".format(gseDir, gse, index)
				eData.to_csv(exprDataFile, sep='\t', index=False, encoding='utf-8')
				metaData.to_csv(metaDataFile, sep='\t', index=True, encoding='utf-8')
				exprDataFiles.append(exprDataFile)
				metaDataFiles.append(metaDataFile)

		# Clean soft.gz files
		subprocess.Popen("rm {0}/*.soft.gz".format(gseDir), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
		# Return metadata and expression files
		return exprDataFiles, metaDataFiles

##--------------------------------------------------------------##
##	Perform bioinformatics analysis								##
##--------------------------------------------------------------##
@retry(stop=stop_after_attempt(2), wait=wait_random(min=1, max=2), retry_error_callback=return_last_value)
def PerformAnalysis(analysis, listParams, logFile):
	# importing analysisCommand var
	global analysisCommand
	# Building the command for the analysis
	KeepLog(logFile, "\t++ Starting {0} \n".format(analysis))
	cmd = analysisCommand[analysis].format(*listParams)
	# KeepLog(logFile, "Starting {0} \n".format(cmd))
	# Launch the command in silent mode
	subprocess.check_call(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=60)
	# KeepLog(logFile, "\t++ Completed {0} \n".format(analysis))

##--------------------------------------------------------------##
##	Download Official MeSH terms								##
##--------------------------------------------------------------##
@retry(stop=stop_after_attempt(2), wait=wait_random(min=1, max=2))
def DownloadOfficialMesh(currentYear, srcFolder, logFile):
	# we first try to download the last/current version of MeSH descriptor from NCBI website

	# It's not necessary to download the file if it's already present
	if (os.path.isfile("{0}/{1}Mesh.xml".format(srcFolder, currentYear))):
		KeepLog(logFile, "\t++ MeSH library for year {0} already downloaded\n".format(currentYear))
	else: # We need to download the file
		try:
			KeepLog(logFile, "\t++ Downloading MeSH library for year {0} \n".format(currentYear))
			downloadCmd = "wget -O {0}/{2}MeSH.xml {1}/desc{2}.xml".format(srcFolder,officialMeshUrl,currentYear)
			subprocess.Popen(downloadCmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
		except:
			KeepLog(logFile, "\t++ {0} not available, downloading MeSH library for previous \n".format(currentYear))
			currentYear = currentYear - 1 # Going back in time of 1 year...
			downloadCmd = "wget -O {0}/{2}MeSH.xml {1}/desc{2}.xml".format(srcFolder,officialMeshUrl,currentYear)
			subprocess.Popen(downloadCmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()

	# Once MeSH library is downloaded, we need to extract informations from it
		# Initialising pandas dataframe for storing meSH info
	meshInfo = pd.DataFrame(columns=["UID","mesh","level","depthLevel"])
		# Loading XML file into a variable
	meshXML = et.parse("{0}/{1}MeSH.xml".format(srcFolder, currentYear))
		# Get root of XML file
	meshXMLRoot = meshXML.getroot()
		# Iterating over MeSH xml and extract informations
	i = 0
	for node in meshXMLRoot:
		try:
			meshUID = node.find("DescriptorUI").text
			meshTerm = node.find("DescriptorName/String").text
			meshLevel = node.find("TreeNumberList/TreeNumber").text
			meshDepthLevel = len(meshLevel.split("."))
			# Update pandas dataframe
			KeepLog(logFile, "\t++ Parsing {0} -- {1}/{2} \n".format(meshUID, i, len(meshXMLRoot)))
			meshInfo.loc[i] = [meshUID, meshTerm, meshLevel, meshDepthLevel]
			i = i+1
		except:
			next

	# Clean downloaded files
	# subprocess.Popen("rm {0}/{1}MeSH.xml".format(srcFolder, currentYear), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	# Return pandas dataframe
	return meshInfo
