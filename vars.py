#! /usr/bin/env python3

##################################################
## SMAC
## File containing the variables for the analysis

## Author: Stefano Pirro, PhD (aka wynstep)
## Version: 2.0
## Maintainer: Stefano Pirro, PhD (aka wynstep)
## Email: stefano.pirro@uniroma2.it, s.pirro@mir-nat.com
## Status: Active
##################################################

##----------------------------------##
##		Folders 					##
##----------------------------------##

rootFolder = "/SMAC2"
scriptsFolder = rootFolder + "/scripts/"
srcFolder = rootFolder + "/src/"
	# Output sub-folders
outRootFolder = "{0}/{1}"
outGraphDir = "{0}/graphs"
outGSEDir = "{0}/GSE"

##----------------------------------##
##	Variables						##
##----------------------------------##
pubmedIDs = [] # List containing the PMIDs that will be retrieved from PubMed
pPMIDs = [] # List containing PMIDs already downloaded and analysed (recovery purposes)
pNuclInfo = [] # List containing Gene mapping details already downloaded and analysed (recovery purposes)
gdsCodes = [] # List containing GDS codes from Gene Expression Omnibus
geneNames = [] # List containing Gene names from PMIDs
geneScores = [] # List containing Genes scores from PMIDs
meshHeadingsContainer = [] # List containing all meSH headings from the analysis
officialMeshUrl = "ftp://nlmpubs.nlm.nih.gov/online/mesh/MESH_FILES/xmlmesh/"
procCont = [] # Container for parallel processing
# Number of elements to get into Entrez query
sizeQuery = { "RetrievePMIDs": 10000, \
			"chunkPMIDs": 1000, \
			"RetrievePMIDInfo": 1000, \
			"RetrieveGeoInfoFromPMID": 1000, \
			"RetrieveMeshInfo": 1000, \
			"PMID2Gene": 1000 }
timoutThr = 120 # Time limit for timout function (seconds)
papersInfo = [] # List containing informations of retrieved papers
geoInfoContainer = [] # List containing informations of GEO datasets
defaultAnalyses = "pca,gene_expression"
analysisCommand = {
	"pca": "Rscript " + scriptsFolder + "/pca.R -e {0} -t {1} -i {2} -n {3} -g {4} -o {5}",
	"receptor_status": "Rscript " + scriptsFolder + "/rStatus.R -e {0} -t {1} -i {2} -n {3} -g {4} -o {5}",
	"tumour_purity": "Rscript " + scriptsFolder + "/tPurity.R -e {0} -t {1} -i {2} -n {3} -g {4} -o {5}",
	"molecular_classification": "Rscript " + scriptsFolder + "/mClass.R -e {0} -t {1} -i {2} -n {3} -g {4} -o {5}",
	"gene_expression": "Rscript " + scriptsFolder + "/gExpr.R -e {0} -t {1} -i {2} -n {3} -g {4} -o {5}",
	"gene_network": "Rscript " + scriptsFolder +"/gNet.R -i {0} -l '{1}' -n {2} -o {3}",
	"rra": "Rscript " + scriptsFolder +"/rra.R -a {0} -b {1} -c {2} -l {3} -n {4} -g {5} -o {6}",
	"functional_enrichment": "Rscript " + scriptsFolder + "/enrichment.R -e {0} -t {1} -i {2} -n {3} -g {4} -o {5}",
	"functional_enrichment_pmid": "Rscript " + scriptsFolder +"/enrichment_pmidGenes.R -i {0} -l '{1}' -n {2} -o {3}"
}

## Files ##
logFile = "{0}/{1}.log"
pmidFile = "{0}/{1}.pmid.txt"
geoInfoFile = "{0}/{1}.geoDS.tsv"
nuclInfoFile = "{0}/{1}.nucleotide.tsv"
nuclPmidFile = "{0}/{1}.nucl_pmid.tsv"
nuclPmidFile_aggregated = "{0}/{1}.nucl_pmid.aggregated.tsv"
geneInfoFile = "{0}/{1}.genes.tsv"
papersInfoFile = "{0}/{1}.papersInfo.tsv"
meshTermsFile = "{0}/{1}.mesh.txt"
meshTermsDetailsFile = "{0}/{1}.mesh.{2}.details.pkl"
meshLAboundanceFile = "{0}/{1}.meshLAboundance.txt"
meshGAboundanceFile = "{0}/{1}.meshGAboundance.txt"
meshLevelFile = "{0}/{1}.meshLevel.txt"
