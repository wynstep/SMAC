# File to create gene list from eData

# importing libraries
import os, fnmatch
import argparse # library for reading arguments

# reading arguments
parser = argparse.ArgumentParser(description='')
parser.add_argument('-d','--dir', help='root directory to search for eData files', required=True)
args = vars(parser.parse_args())
root_dir= args["dir"]

## Custom function to list all the files in a root directory, given a pattern
def recursive_glob(treeroot, pattern):
    results = []
    for base, dirs, files in os.walk(treeroot):
        goodfiles = fnmatch.filter(files, pattern)
        results.extend(os.path.join(base, f) for f in goodfiles)
    return results

all_eData_files = recursive_glob(root_dir, "eData.*.tsv")

for eData in all_eData_files:
	eData_basename = eData.split("/")[-1]
	eData_directory = eData.replace(eData_basename,'')
	eData_num = eData_basename.split(".")[1]

	# apply system command to create gene list
	command = "cat {0} | cut -f1 > {1}gene.{2}.list".format(eData,eData_directory,eData_num)
	os.system(command)

