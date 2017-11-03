#!/usr/bin/env python

import sys
import math
import numpy as np
import shapeLibrary as sl
import argparse

nucl = ("A","C","G","T")

##################  VARIOUS FUNCTIONS  #################

##################  MAIN FUNCTION  #################

def main():

	#Creating parser
	parser = argparse.ArgumentParser(description='Reads a list of sequences and computes the mean shape profile using a k-mer table. Outputs one position in the profile per line. The k-mer table can have multiple columns, leading to multiple columns of output.')
	parser.add_argument('seqFile', metavar='seq.lst', help='Text file containing one sequence per line ("-" gives STDIN)')
	parser.add_argument('kmerFile', metavar='kmerValue.csv', help='Comma-separated kMer file. The first column contains k-mers, the following columns contain the associated value(s).')
	parser.add_argument("--scoreN", help="Treats N as the average of A, C, G, and T.", action="store_true")
	parser.add_argument("--header", help="First line in k-mer file is header.", action="store_true")
#	parser.add_argument("--verbose", help="Increase output verbosity", action="store_true")
	args = parser.parse_args()

	#Parses kmer file
	(kMers, values, k, nCol, header) = sl.readKMerTable(args.kmerFile, args.header)
	kMerTable 			= dict([(kMers[i], values[i]) for i in range(len(kMers))])

	if args.scoreN: #Adds wildcard "N" characters to the k-mer table by averaging over "A", "C, "G", "T"
		for x in range(k): #Successively adds "N" to each position (in combination with previously added Ns)
			currentKeys = kMerTable.keys()
			for key in currentKeys:
				newKey 	= key[:x]+"N"+key[x+1:]
				if newKey in kMerTable:
					continue
				else:
					kMerTable[newKey] = sum([ kMerTable[key[:x]+n+key[x+1:]] for n in nucl ]) / len(nucl)
		

	#Computes the mean value
	nSeq				= 0
	seqSum				= None
	L					= None

	#Determines where to read the sequence file
	if args.seqFile == "-":
		f = sys.stdin
	else:
		f = open(args.seqFile)

	#Loops over sequences
	for l in f:
		#Makes sure the sequences have equal length (and sets up seqSum the first round0
		if L is None:
			L		= len(l.rstrip())
			seqSum	= [np.zeros(nCol) for i in range(L-k+1)]
		elif L != len(l.rstrip()):
			sl.err('All sequences must be of equal length.')

		nSeq		+=1
		seq			= l.rstrip()

		for i in range(L-k+1):
			seqSum[i] +=  kMerTable[l[i:i+k]]

	#Prints the mean profile
	if args.header:
		print ",".join(header)
	for i in range(L-k+1):
		print ",".join([ "%f"%di for di in (seqSum[i]/nSeq)])
		
##################  END #################

main()
