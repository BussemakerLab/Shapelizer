#!/usr/bin/env python

import sys
import math
import numpy as np
import scipy.optimize as op
import shapeLibrary as sl
import argparse
from string import maketrans

nucl = ("A","C","G","T")

##################  VARIOUS FUNCTIONS  #################

def scoreKMer(seq, m):
	"""Scores a k-mer using a tensor containing a scoring matrix and dinucleotide interactions."""
	s = 0
	numSeq = [nucl.index(b) for b in seq]
	#Mononucleotide contribution
	s += sum([ m[0][numSeq[i]][i] for i in range(len(seq)) ])
	#Dinucleotide contribution
	s += sum([ m[1][4*numSeq[i]+numSeq[i+1]][i] for i in range(len(seq)-1) ])
	
	return s
		
##################  MAIN FUNCTION  #################

def main():

	#Creating parser
	parser = argparse.ArgumentParser(description='Generates a random k-mer table from an input table using either the "matched complexity" method (default) or by permuting the input table.')
	parser.add_argument('kmerFile', metavar='kmerTable.csv', help='Comma-separated kMer table. (COL 1) = kmer, (COL 2) = values')
	parser.add_argument("-s", metavar='seed', help="Seed for numpy.random", default=None)
	parser.add_argument("--verbose", help="Increase output verbosity", action="store_true")
	parser.add_argument("--symmetric", help="k-mer table is reverse-complement symmetric.", action="store_true")
	parser.add_argument("--header", help="First line in kmer file is header.", action="store_true")
	parser.add_argument("--permute", help="Permuts the input table", action="store_true")
	args = parser.parse_args()

	#Sets seed.
	if args.s is not None:
		np.random.seed(int(args.s))

	#Reads k-mer table
	(kMers, values, k, nCol, header)		= sl.readKMerTable(args.kmerFile, args.header)
	if len(values[0]) > 1:
		sl.err("Current implementation can only generate a single random column")

	if args.permute: #Permutes k-mer table
		if args.symmetric:
			#Identifies reverse-complement pairs of sequences that shouldbe held out 
			kmerPairs			= {}
			symValues			= {}
			trantab				= maketrans("ACGT", "TGCA")
			for i in range(len(kMers)):
				km				= kMers[i]
				rcKm				= km[::-1].translate(trantab)
				if rcKm not in kmerPairs:
					kmerPairs[kMers[i]]	= (i, kMers.index(rcKm))
					symValues[kMers[i]]     = (values[i] + values[kMers.index(rcKm)]) / 2

			#Creates list of inital and permuted kmers
			inKmers                         = kmerPairs.keys()
			pKmers                          = np.random.permutation(inKmers)
			newValues                       = [0.] * len(kMers) 
			for i in range(len(inKmers)):
				km                      = inKmers[i]
				rcKm                    = km[::-1].translate(trantab)
				permutedValue           = symValues[pKmers[i]]

				newValues[kMers.index(km)]      = permutedValue
				newValues[kMers.index(rcKm)]    = permutedValue 

		else:
			newValues						=  np.random.permutation(values[:,0])
	else: #Generates random mono+di model with matched conditional variance

		# 1. Computes the expected conditional variance of 'true' k-mer table
		conditionalVariance 				= np.array([ [ np.mean([ np.var([ values[j][0] for j in range(len(kMers)) if kMers[j][i1] == nucl[n1] and kMers[j][i2] == nucl[n2] ]) for n1 in range(4) for n2 in range(4) if  not (i1==i2 and n1!=n2)]) for i2 in range(k)] for i1 in range(k) ])

		# 2. Create a design matrix used to generate random k-mer table
		# 2.1 Creates (independent) mononucleotide matrices:
		sigList								= []
		matrixList							= []
		for i in range(k):
			#Generates new mononucleotide-matrix using uniform random numbers
			rv								= np.random.rand(4)
			newMonoMatrix					= np.array([ [0.]*(i)+[rv[y]-np.mean(rv)]+[0.]*(k-i-1) for y in range(4)])
			#Generates a "signature" vector used to make sure we don't add the same degrees of freedom twice (for symmetric matrices)
			newMonoSignature				= np.array([0]*(i)+[1]+[0]*(k-i-1))
	
			#Symmerizes matrix and signature if appropriate.
			if args.symmetric:
				newMonoMatrix				= (newMonoMatrix +  newMonoMatrix[::-1,::-1]) / 2
				newMonoSignature			=  newMonoSignature + newMonoSignature[::-1]

			#Checks if the current signature allready has been added.
			sig = "mono"+"".join(["%d"%si for si in newMonoSignature])
		
			#Saves the new matrix if an equivalent matrix has not been saved before
			if sig in sigList:
				continue
			else:
				sigList						+= [sig]
				matrixList					+= [ (newMonoMatrix, np.zeros((16,k-1)))]
	

		# 2.2 Creates (independent) dinucleotide matrices:
		#matrix used to reverse-complement a 16-entry dinucleotide vector
		rcDiMatrix					= np.array([ [ int(j==4*(3-n2)+(3-n1)) for n1 in range(4) for n2 in range(4) ] for j in range(16) ])
		for i in range(k-1):
			#Generates new mononucleotide-matrix using uniform random numbers
			rv								= np.random.rand(16)
			newDiMatrix						= np.array([ [0.]*(i)+[rv[y]-np.mean(rv)]+[0.]*(k-i-2) for y in range(16)])
			#Generates a "signature" vector used to make sure we don't add the same degrees of freedom twice
			newDiSignature					= np.array([0]*(i)+[1]+[0]*(k-i-2))

			#Symmerizes matrix and signature if appropriate.
			if args.symmetric:

				newDiMatrix					= (newDiMatrix    + rcDiMatrix.dot(newDiMatrix[:,::-1])) / 2
				newDiSignature				=  newDiSignature + newDiSignature[::-1]

			#Checks if the current signature allready has been added.
			sig = "di"+"".join(["%d"%si for si in newDiSignature])
		
			#Saves the new matrix if an equivalent matrix has not been saved before
			if sig in sigList:
				continue
			else:
				sigList						+= [sig]
				matrixList					+= [ (np.zeros((4,k)), newDiMatrix) ]

		# 2.3 Computes design matrix by scoring each k-mer using each independent mono/di-nucleotide matrix 
		X 									= np.array([ [ scoreKMer(kMer, m) for m in matrixList ]for kMer in kMers])	

		# 3. Computes the conditional covariance between pairs of matrix models
		d 									= [[ np.array([[np.mean([ np.cov(np.array([ [X[j][a], X[j][b]] 
			for j in range(len(kMers)) if kMers[j][i1] == nucl[n1] and kMers[j][i2] == nucl[n2] ]).transpose())[0,1] 
				for n1 in range(4) for n2 in range(4) if  not (i1==i2 and n1!=n2)]) 
					for a in range(len(matrixList))] for b in range(len(matrixList)) ])
						for i1 in range(k)] for i2 in range(k) ]

		# 4. Finds a combination of matrix models that minimizes the L2-error in the expected conditional variance.
		# 4.1 Loss function
		f  									= lambda v: np.sum([   (v.dot(d[i1][i2]).dot(v)-conditionalVariance[i1,i2])**2               for i1 in range(k) for i2 in range(k) ])
		# 4.2 Gradient of loss function
		df 									= lambda v: np.sum([ 4*(v.dot(d[i1][i2]).dot(v)-conditionalVariance[i1,i2])*v.dot(d[i1][i2]) for i1 in range(k) for i2 in range(k) ], axis=0)
		# 4.3 Initial seed
		x0 									= np.random.rand(len(matrixList))
		# 4.4 Minimizes the loss funciton using L-BFGS
		res 								= op.minimize(f, x0, args=(), method='L-BFGS-B',    jac=df, options={'disp': False, 'maxiter':1000})
		# 4.5 Computes new values
		newValues 							= np.array([Xi.dot(res.x) for Xi in X])

		if args.verbose:
			# Writes conditional variance matrices
			sl.disp("Conditional variance matrix:")
			sl.printMatrix(conditionalVariance, sys.stderr)

			sl.disp("Conditional variance in random model:")
			newConditionalVariance 			= np.array([ [ np.mean([ np.var([ newValues[j] for j in range(len(kMers)) if kMers[j][i1] == nucl[n1] and kMers[j][i2] == nucl[n2] ]) for n1 in range(4) for n2 in range(4) if  not (i1==i2 and n1!=n2)]) for i2 in range(k)] for i1 in range(k) ])
			sl.printMatrix(newConditionalVariance, sys.stderr)

	#Writes the random matrix to STDOUT
	for i in range(len(newValues)):
		print "%s,%f"%(kMers[i], newValues[i])
	

##################  END #################

main()
