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
def kMerToVector(kMer, mono, di, all):
	"""Converts a kMer to a line in a design matrix"""

	
	if mono or all:
		#Computes mononucleotide representation
		monoRep			= dict([ (nucl[i], np.identity(4)[i].tolist()) for i in range(len(nucl))] )
		x				= []
		for k in kMer:
			x			+=  monoRep[k]

		if all:
			#Takes outer product to get all-by-all representation
			x= np.tensordot(np.array(x),np.array(x), axes=0).flatten().tolist()

	else:
		#Computes dinucleotide representation
		diNucl				= [ n1+n2 for n1 in nucl for n2 in nucl ]
		diRep				= dict([ (diNucl[i], np.identity(16)[i].tolist()) for i in range(len(diNucl))] )
		x					= []
		for i in range(len(kMer)-1):
			x				+= diRep[kMer[i:i+2]]

	return x
	
##################  MAIN FUNCTION  #################

def main():

	#Creating parser
	parser = argparse.ArgumentParser(description='Reads a k-mer model and evaluates in for all kmers.')
	parser.add_argument('modelFile', metavar='motel.tsv', help='Model file')
	parser.add_argument("--verbose", help="Increase output verbosity", action="store_true")
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
