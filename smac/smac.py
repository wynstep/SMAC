#!/usr/bin/env python
#title           :smac.py
#description     :This will execute the main SMAC script
#author          :s.pirro
#date            :20180123
#version         :0.1
#usage           :python smac.py
#notes           :
#==============================================================================

# importing libraries
import os, sys
import argparse # library for reading arguments
from vars import * # importing custom variables

# reading arguments
parser = argparse.ArgumentParser(description='Welcome to SMAC')
parser.add_argument('-t','--terms', help='terms for performing the analysis', required=False, nargs='+')
parser.add_argument('-p','--pmids', help='file with list of pmids to download and analyse', required=False)
parser.add_argument('-l','--limit', help='max number of publications to retrieve', required=False)
parser.add_argument('-a','--analysis', help='analysis to perform (separated by comma) [pca,receptor_status,molecular_classification,tumour_purity,gene_expression]', required=False)
parser.add_argument('-s','--skip', help='do you want to skip GSE download/analysis [1=yes, skip the analysis, 0=no]', choices=['0','1'], required=False)
parser.add_argument('-e','--email', help='email for retrieving the papers', required=True)

#==============================================================================
#	Checking arguments
#==============================================================================

args = vars(parser.parse_args())
output_folder = "results"

# Checking if the user wants to skip the GSE download/analysis
if args["skip"] is None:
	args["skip"] = 1

# In case the user does not want to skip, check if the analysis have been selected
if args["skip"] == 0:
	if args["analysis"] is None:
		sys.exit(exec_2_err)
	else:
		submitted_analysis = args["analysis"].rstrip().split(",")
		matched_analysis = list(set(valid_analysis) & set(submitted_analysis))
else:
	# if the skip flag is active, we do not care about the analysis, we just set a default value
	matched_analysis = ["pca"]


#==============================================================================
#	Running commands
#==============================================================================

# Here we intersect the submitted analyses with the available and valid. The script will continue just if the user submitted valid analyses, otherwise it will stop
if (len(matched_analysis) > 0):
	# Checking presence of one parameter between -t and -p
	if (args["terms"] is None and args["pmids"] is not None): # Here we check the pmid list
		# We call the subscript specific for the pmids
		os.system("python3 smac_pmid.py -p {0} -a {1} -e {2} -o {3} -s {4}".format(args["pmids"],",".join(matched_analysis),args["email"],output_folder,args["skip"]))
	elif (args["pmids"] is None and args["terms"] is not None): # Here we check the terms
		# We call the subscript specific for the terms
		os.system("python3 smac_term.py -t {0} -l {1} -a {2} -e {3} -o {4} -s {5}".format(" ".join(args["terms"]),args["limit"],",".join(matched_analysis),args["email"],output_folder,args["skip"]))
	else:
		# Here both terms and pmids were not set
		print(exec_1_err)
else:
	# there is no valid analysis to submit
	print(exec_2_err)
