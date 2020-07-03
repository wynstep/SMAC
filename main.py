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
import os, sys, argparse, subprocess
from vars import *
from functions import *

	## Reading arguments for the analysis
parser = argparse.ArgumentParser()
parser.add_argument('-t','--type', help='List type of SMAC analysis to perform <term, pmid, gene>', required=True)
parser.add_argument('-p','--params', help='Arguments associated to the specific type', required=True, nargs=argparse.REMAINDER)
args = parser.parse_args()

# Define type of analysis
typeAnalysis = args.type
params = " ".join(args.params)

# Clearing terminal and print welcome message
os.environ['TERM'] = 'xterm'
os.system("clear")
initialMessage = "######### SMAC 2.0 #########\n#\n# Type analysis: {0}\n# Parameters: {1}\n#\n############################\n\n".format(typeAnalysis, params)
terminalPrint(initialMessage)

# Launch a different script according to the type of analysis
if (typeAnalysis == "gene"):
	cmd = "python3 main-gn.py {0}".format(params)
elif (typeAnalysis == "term"):
	cmd = "python3 main-term.py {0}".format(params)
elif (typeAnalysis == "pmid"):
	cmd = "python3 main-pmid.py {0}".format(params)
else:
	print("[ERROR] You need to declare a value between <gene, term, pmid>")
	sys.exit()

# Launching analysis according to the type
subprocess.Popen(cmd, shell=True).communicate()
