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
	parser = argparse.ArgumentParser(description='Uses the kmer table to score sequence and computes the mean value. One value per line.')
	parser.add_argument('seqFile', metavar='seq.lst', help='Text file containing one sequence per line ("-" gives STDIN)')
	parser.add_argument('kmerFile', metavar='kmerValue.csv', help='kMer file. (COL 1) = kmer, (COL 2) = values')
	parser.add_argument("--scoreN", help="Treats N as the average of A, C, G, and T.", action="store_true")
	parser.add_argument("--header", help="First line in kmer file is header.", action="store_true")
	parser.add_argument("--verbose", help="Increase output verbosity", action="store_true")
	args = parser.parse_args()

	#Parses kmer file
	kMerTable			= {}
	header				= None
	k					= None
	nCol				= None
	with open(args.kmerFile) as f:
		if args.header:
			header		= f.readline().rstrip().split(',')[1:]
			
		for l in f:
			d			= l.rstrip().split(",")

			if k is None:
				k		= len(d[0])
				nCol	= len(d)-1
			elif k != len(d[0]):
				sf.err('All kmers must be of equal length')

			try:
				kMerTable[d[0]] = np.array([ float(di) for di in d[1:] ])
			except ValueError:
				sl.err("Table entries must be float numbers. LINE='%s'."%(l.rstrip()))

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
