#! /usr/bin/env python3

################################################################################
#
#   File name: createWordCloud.py
#
#   Authors: Stefano Pirro', PhD ( s.pirro@mir-nat.com)
#
#
#		Description: Script for creating wordcloud MeSH terms
#
################################################################################

	# Importing libraries, vars and functions
import os, sys, argparse
import pandas as pd
import numpy as np
import matplotlib as mpl
if os.environ.get('DISPLAY','') == '':
	print('no display found. Using non-interactive Agg backend')
	mpl.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
from wordcloud import WordCloud, STOPWORDS, ImageColorGenerator

	## Reading arguments for the analysis
parser = argparse.ArgumentParser()
parser.add_argument("-a", "--statFile", action="store", required=True, help="TSV containing meSH details")
parser.add_argument("-l", "--label", action="store", required=True, help="label for plot")
parser.add_argument("-o", "--outputFile", action="store", required=True, help="Wordcloud file name")
args = parser.parse_args()
## Create word bag
wordBag = {}

## Loading meshFile
statData = pd.read_csv(args.statFile, sep="\t")

if args.label == "mirnas" or args.label == "gene":
	aggregatedData = statData.groupby(["geneName"]).agg(['count'])
	for gene, row in aggregatedData.iterrows():
		wordBag[gene] = row["PMID"]["count"]
else:
	for index, row in statData.iterrows():
		wordBag[row[0]] = row[1]

	## Create image mask
imageMaskFile = "/media/extData/MAP/html/images/wcMask.png"
imageMask = np.array(Image.open(imageMaskFile))

	# Generate wordcloud
wordcloud = WordCloud(max_font_size=90, max_words=1000, random_state=42, background_color="white", mask=imageMask).generate_from_frequencies(frequencies=wordBag)
# create coloring from image
imageColors = ImageColorGenerator(imageMask)
plt.figure(figsize=[7,7])
plt.imshow(wordcloud.recolor(color_func=imageColors), interpolation="bilinear")
plt.axis("off")
plt.show()
wordcloud.to_file(args.outputFile)
