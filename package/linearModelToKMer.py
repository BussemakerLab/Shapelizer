#!/usr/bin/env python

import sys
import math
import numpy as np
import scipy as sp
import scipy.stats as st
import shapeLibrary as sl
import argparse
from string import maketrans

nucl = ("A","C","G","T")

##################  VARIOUS FUNCTIONS  #################
	
##################  MAIN FUNCTION  #################

def main():

	#Creating parser
	parser = argparse.ArgumentParser(description='Reads a k-mer model and evaluates in for all kmers. Outputs CSV file.')
	parser.add_argument('modelFile', metavar='motel.tsv', help='Model file. Same format as output from "kMerLinearRegression.py --betas"')
#	parser.add_argument("--verbose", help="Increase output verbosity", action="store_true")
	args = parser.parse_args()

	model        = sl.loadScoringMatrix(args.modelFile, readIntercept=True)

	k = len(model[0])

	#Generates a list of kMers
	kMers = [""]
	for x in range(k):
		temp     = []
		for ki in kMers:
			temp += [ ki+n for n in nucl ]
		kMers = temp

	for ki in kMers:
		print "%s,%f"%(ki, sl.scoreSeq(ki, model))


##################  END #################

main()
