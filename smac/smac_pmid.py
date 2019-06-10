#!/usr/bin/env python
#title           :smac.py
#description     :This will execute the SMAC script for pmids
#author          :s.pirro
#date            :20180123
#version         :0.1
#usage           :python smac_pmid.py -t {pmids} -l {limit_pubs} -a {analyses} -e {email} -o {output_dir}
#notes           :
#==============================================================================

# importing libraries
import os, sys
import random, string # library for generating random code
import argparse # library for reading arguments
from functions import * # importing custom functions
from vars import * # importing custom variables

# reading arguments
parser = argparse.ArgumentParser(description='Welcome to SMAC')
parser.add_argument('-p','--pmids', help='file with list of pmids to download and analyse', required=True)
parser.add_argument('-a','--analysis', help='analysis to perform (separated by comma) [pca,receptor_status,molecular_classification,tumour_purity,gene_expression]', required=True)
parser.add_argument('-n','--ngenes', help='Number of top differentially-expressed genes to plot in the gene expression analysis', required=False)
parser.add_argument('-e','--email', help='email for retrieving the papers', required=True)
parser.add_argument('-o','--output', help='output for results', required=True)
parser.add_argument('-s','--skip', type=int, help='do you want to skip GSE download/analysis [1=yes, skip the analysis, 0=no]', choices=[0,1], required=True)
args = vars(parser.parse_args())

# initialising the output folder
output_folder = args["output"]

### creating directories where to store the results
# Initialising a random hexcode that will be the unique id for the analysis
random_code = ''.join(random.choices(string.ascii_lowercase + string.digits, k=6))

# create subfolders (if they not exist)
root_folder = "{0}/{1}".format(output_folder, random_code)
data_folder = "{0}/data".format(root_folder)
geo_folder = "{0}/geo_datasets".format(data_folder)
graph_folder = "{0}/graphs".format(root_folder)
ent_min_folder = "{0}/entropy_minimisation".format(graph_folder)
if not os.path.exists(root_folder):
	os.makedirs(root_folder)
	os.makedirs(data_folder)
	os.makedirs(geo_folder)
	os.makedirs(graph_folder)
	os.makedirs(ent_min_folder)

### Starting the execution of the script
os.system("clear")
# print initial message
print(welcome_pmid % (args["pmids"], args["email"], random_code))

# collecting the pmid items into a list
paper_list = []
with open(args["pmids"], "r") as pmid_io_stream:
	for line in pmid_io_stream:
		paper_list.append(line.rstrip())

# Retrieve info and MeSH Headers from papers
all_infos = RetrieveInfoFromPMID(paper_list, data_folder, args["email"])
mesh_terms_no_unique = all_infos[1]
mesh_terms = all_infos[2]
papers_info = all_infos[3]

# Retrieve GEO datasets from papers and analysing -- if skip=0 is selected
if (args["skip"] == 0):
	for pmid in papers_info.keys():
		# total number of gses associated to the pmid
		tot_gse = len(papers_info[pmid]["gse_codes"])
		# create pmid_folder (just for pmids with GEO data associated)
		if tot_gse>0:
			pmid_folder = "{0}/{1}".format(geo_folder, pmid)
			# check if folder already exists
			if not os.path.exists(pmid_folder):
				os.makedirs(pmid_folder)

			for k,gse in enumerate(papers_info[pmid]["gse_codes"]):
				# creating gse folder (subfolder of pmid)
				gse_folder = "{0}/GSE{1}".format(pmid_folder, gse)
				# check if folder already exists
				if not os.path.exists(gse_folder):
					os.makedirs(gse_folder)

				# Download GEO dataset
				PrintInline("	+-- Downloading GSE{0} -- PMID: {1} -- {2}/{3}".format(gse, pmid, k+1, tot_gse))
				DownloadGEODataset(gse, gse_folder)
				## Create target file for each pData downloaded
				CreateTargetFile(gse_folder)

				# listing the target files in the gse folder
				target_files = glob.glob("{0}/target.*.tsv".format(gse_folder))

				# Perform the selected analyses for each expression data and target file
				print("\n## Performing selected analyses...")
				# Iterating into target files and retrieve the right number of the analysis
				for t_File in target_files:
					cont = re.search('.+([0-9]).tsv', t_File).group(1)
					expression_file = "{0}/eData.{1}.tsv".format(gse_folder,str(cont))
					PerformAnalyses(gse_folder, t_File, expression_file, args["analysis"], args["ngenes"], str(cont))

				# Update papers infos with the performed analyses
				print("\n## Updating infos (with performed analyses)...")
				papers_info = UpdatePapersInfo(papers_info, pmid, gse, gse_folder, args["analysis"], data_folder)

				### Update 2019-05-21 -- Before continuing we want to be sure that the folder we created is not empty, otherwise the folder will be removed
				CheckEmptyFolder(gse_folder)

			### Update 2019-05-22 -- Before continuing we want to be sure that the folder we created is not empty, otherwise the folder will be removed
			CheckEmptyFolder(pmid_folder)

# Save final publications
print("## Save literature report")
CreateLiteratureReport(papers_info, data_folder)

# Retrieve genes from papers
print("## Retrieving genes related to the publications...")
gene_list = RetrieveGenesFromPMID(paper_list, data_folder)
# Map genes on Mentha interactome database and create gene networks
if (len(gene_list) > 0):
	print("## Create gene network for selected genes...")
	gene_list_file = "{0}/gene_list.txt".format(data_folder)
	network_file = "{0}/mentha.txt".format(src_folder)
	command = "Rscript scripts/PubMed/CreateNetwork.R --list {0} --output {1} --network {2}".format(gene_list_file, root_folder, network_file)
	os.system(command)
else:
	print(exec_5_err % "Network Creation", end="\r")

# Sort lists of MeSH terms
if (len(mesh_terms) > 0):
	print("## Sorting MeSH terms...")
	sorted_by_level = SortMesh("level", mesh_terms, data_folder)
	sorted_by_local_aboundance = SortMesh("laboundance", mesh_terms_no_unique, data_folder)
	sorted_by_global_aboundance = SortMesh("gaboundance", mesh_terms, data_folder)
	# Combine lists of MeSH terms
	print("## Aggregating MeSH terms...")
	command = "Rscript scripts/PubMed/AggregateLists.R --lists {0},{1},{2} --output {3}".format(sorted_by_level,sorted_by_local_aboundance,sorted_by_global_aboundance,root_folder)
	os.system(command)
    # create animated gif from generated images
	command = "convert -delay 20 -loop 0 $(find {0} -name 'update_*.png' -print0 | sort -zV | xargs -r0 echo) {1}/entropy_minimisation.gif".format(ent_min_folder,graph_folder)
	os.system(command)
	# Calculate p-value fluctuations
	print("## Calculating p-value and kendall-tau fluctuations")
	command = "Rscript scripts/PubMed/CalculatePvalue.R --list {0}/aggregated_mesh.txt --output {1}".format(data_folder, root_folder)
	os.system(command)
else:
	print(exec_5_err % "MeSH aggregation")
