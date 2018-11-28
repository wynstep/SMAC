#!/usr/bin/env python
#title           :vars.py
#description     :This file contains all the variables for smac
#author          :s.pirro
#date            :20180123
#version         :0.1
#==============================================================================

## Error messages (we name by code)
exec_1_err = """Sorry, I cannot exec the SMAC program because you did not select any term or pmid file.
              Please run 'python smac.py -h' for details"""
exec_2_err = """Sorry, I cannot exec the SMAC program because you did not submit any valid analysis.
              Please run 'python smac.py -h' for details"""
exec_3_err = """Error -- GSE%s folder will be deleted. For possible reasons, please refer to the documentation"""
exec_4_err = """ERROR -- %s is not a valid analysis. Please refer to the documentation for details"""
exec_5_err = """ERROR in performing %s -- The analysis will continue"""

## General messages
welcome_pmid = """#######################################################
#  SMAC -- The Smart Automatic Classification Method  #
#######################################################
Parameter: %s
email: %s
unique code for the analysis: %s

"""

## List of valid analyses the user can perform
valid_analysis = ["tumour_purity","molecular_classification","receptor_status","pca","gene_expression"]

## List of background terms to remove from samples descriptions
background_terms = ["samples?\s\w+.*?(?=\s)", "replicates?\s\w+.*?(?=\s)"]

## Folders and file names (fixed)
src_folder = "src"

