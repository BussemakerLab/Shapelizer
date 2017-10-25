#!/usr/bin/env python

import math
import numpy as np
import shapeLibrary as sl
import argparse

nucl = ("A","C","G","T")

##################  VARIOUS FUNCTIONS  #################

def sampleSequence(scoringMatrix, beta, nOut):
	"""Samples sequences using monteCarlo sampling""" 
	k 						= len(scoringMatrix) # Length of scoring matrix
	seq						= [ np.random.randint(0,4) for i in range(k) ] # Current sequence
	E						= sum([ scoringMatrix[i,seq[i]]for i in range(k)]) # Current energy

	seqNr					= 0	
	testNr					= 0
	testList				= None
	output					= [None] * nOut
	outputNr				= 0
	while outputNr < nOut:
	
		#Finds a random (position, base) combination to test
		if testNr%(4*k) == 0 or testList is None:
			testList		= np.random.permutation([ (i,n) for i in range(k) for n in range(4)])
		(iTest, nTest)		= testList[testNr]
		testNr				= (testNr+1)%(4*k)

		if seq[iTest] != nTest:
			seqNr 			+= 1
			deltaE			= scoringMatrix[iTest,nTest] - scoringMatrix[iTest,seq[iTest]]

			#Metropolis-Hastings sampling:
			if beta*deltaE > 0 or np.random.rand() < np.exp(beta*deltaE):
				seq[iTest]	= nTest
				E			+= deltaE
				
				#Saves every 2*k sequence after 10*k burn-in edits
				if seqNr>10*k and seqNr%(2*k) == 0:
				
					output[outputNr]	= (seq[:], E)
					outputNr 			+=1

					
			#Recomputes E every 1000 iterations (to suppress numerical errors, can be removed)
			if seqNr%1000==0:
				E						= sum([ scoringMatrix[i,seq[i]]for i in range(k)]) 

	return output

##################  MAIN FUNCTION  #################

def main():

	#Creating parser
	tDEFAULT			= 13
	eDEFAULT			= 10
	nDEFAULT			= 1000
	parser = argparse.ArgumentParser(description='Samples sequences across DDG bins using Metropolis-Hastings sampling.')
	parser.add_argument('matrixFile', metavar='scoringMatrix.tsv', help='Mononucleotide scoring matrix. (COL 1) = "A,C,G,T", (COL 2..) = -ddG values')
	parser.add_argument('-t', metavar='nTemp', help='Number of temperatures (DEFAULT = %d)'%tDEFAULT, type=int, default=tDEFAULT)
	parser.add_argument('-e', metavar='nEnergyBins', help='Number of energy bins (DEFAULT = %d)'%eDEFAULT, type=int, default=eDEFAULT)
	parser.add_argument('-n', metavar='nSample', help='Number of sampled sequences (DEFAULT = %d)'%nDEFAULT, type=int, default=nDEFAULT)
	parser.add_argument("--verbose", help="Increase output verbosity", action="store_true")
	args = parser.parse_args()

	#1 Reads and processes scoring matrix
	#1.1 - Reads matrix file 
	rawScoringMatrix	= sl.loadScoringMatrix(args.matrixFile)
	if rawScoringMatrix.keys() != [0]:
		sl.err("Matrix file should only contain mononucleotides")

	#1.2. Normalizes so score is in the interval [-1,0]
	scoringMatrix		= np.array([ i-max(i) for i in rawScoringMatrix[0] ]) / sum([ max(i)-min(i) for i in rawScoringMatrix[0] ])
	k            		= len(scoringMatrix)
	dE					= 1. / args.e

	#Samples sequences at each temperature
	betaList 			= 2 * np.log(math.pow(4,k)) *  np.linspace(-1,1,args.t)
	binnedSeq			= [ [ [] for ei in range(args.e) ] for ti in range(args.t) ] 

	for betaI in range(len(betaList)):
		#Samples sequences with using Metropolis-Hastings sampling from distribution p[s] ~ e^(beta*E)
		beta			= betaList[betaI]
		seqList			= sampleSequence( scoringMatrix, beta, 3*args.n )

		for s in seqList:
			#Identifies energy bin
			energyI 	= min(int(np.floor((s[1]+1)*args.e)), args.e-1)
			#In each bin, saves the sequence with probability proportional to e^(-beta*E)
			if (beta>0 and np.random.rand()<np.exp(-beta*( s[1]-(-1+energyI*dE)))) or (beta<=0 and np.random.rand()<np.exp(-beta*(s[1]-(-1+(energyI+1)*dE)))):
					binnedSeq[betaI][energyI] += [s]


	#Builds sequences
	outSeqs=[]
	for iE in range(args.e):
		outSeqs+=[[]]
		for Ti in range(args.t):
			for s in binnedSeq[Ti][iE]:
				outSeqs[iE]+=["".join(nucl[n] for n in s[0])]
		outSeqs[iE] =  np.random.permutation(outSeqs[iE])

	for iE in range(args.e):
		for iS in range(min(args.n,len(outSeqs[iE]))):
			print ">iE=%d,iSeq=%d"%(iE+1,iS+1)
			print outSeqs[iE][iS]
			
#		print ",".join([ "".join(nucl[n] for n in s[0]) for Ti in range(args.t) for s in binnedSeq[Ti][Ei] ])

#DEBUG: Prints all the energies
#	for i in range(len(binnedSeq)):
#		print "\t".join([ ",".join(["%f"%binnedSeq[i][j][x][1] for x in range(len(binnedSeq[i][j]))]) for j in range(len(binnedSeq[i]))])

##################  END #################

main()
