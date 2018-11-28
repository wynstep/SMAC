#!/usr/bin/env python
#title           :BuildMeshDict.py
#description     :This will build MeSH dictionary
#author          :s.pirro
#date            :20180123
#version         :0.1
#usage           :python BuildMeshDict.py --res [mesh_res]
#notes           :
#==============================================================================

# importing libraries
import argparse # library for reading arguments

# reading arguments
parser = argparse.ArgumentParser(description='Welcome to SMAC')
parser.add_argument('-r','--res', help='mesh dictionary to load', required=True)
parser.add_argument('-o','--output', help='output for results', required=True)
args = vars(parser.parse_args())
output_folder = args["output"]

# initialising vars
mesh_dict = {}

# read dictionary file
file_io = open(args["res"], "r").read()
# splitting by fixed string
descriptors = file_io.split("*NEWRECORD")
# iterating and search for Heading and Level
for descriptor in descriptors:
    all_lines = descriptor.split("\n")
    mh = level = ""
    for line in all_lines:
        if "MH = " in line:
            mh = line.split(" = ")[-1]
        elif "MN = " in line:
            level = len(line.split(" = ")[-1].split("."))
    mesh_dict[mh] = level

# writing dictionary into a file
file_io = open("{0}/mesh_dict.txt".format(output_folder), "w")
for mh, level in mesh_dict.iteritems():
    file_io.write("{0}\t{1}\n".format(mh, level))
file_io.close
