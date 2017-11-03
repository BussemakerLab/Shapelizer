#!/usr/bin/env python

import math
import numpy as np
import shapeLibrary as sl
import argparse

nucl = ("A","C","G","T")

##################  VARIOUS FUNCTIONS  #################

def sampleSequence(scoringMatrix, beta, nOut):
	"""Samples nOut sequences using monteCarlo sampling with inverse temperature beta and energy computed using the scoring matrix""" 
	k 						= len(scoringMatrix)                               # Width of scoring matrix
	seq						= [ np.random.randint(0,4) for i in range(k) ]     # Current sequence in numeric representation
	E						= sum([ scoringMatrix[i,seq[i]]for i in range(k)]) # Current energy

	seqNr					= 0
	testNr					= 0
	testList				= None
	output					= [None] * nOut
	outputNr				= 0

	while outputNr < nOut:
		#Generates a randomly permuted list of all possible edits (position, base)  to try. Then tries one at a time.
		if testNr%(4*k) == 0 or testList is None:
			testList		= np.random.permutation([ (i,n) for i in range(k) for n in range(4)])
		#Pulls the next edit from the list of permted edits
		(iTest, nTest)		= testList[testNr]
		testNr				= (testNr+1)%(4*k)

		#If the edit is non-trivial...
		if seq[iTest] != nTest:
			#... computes the chang in energy.
			deltaE			= scoringMatrix[iTest,nTest] - scoringMatrix[iTest,seq[iTest]]
			seqNr 			+= 1

			#Metropolis-Hastings sampling:
			if beta*deltaE > 0 or np.random.rand() < np.exp(beta*deltaE):
				seq[iTest]	= nTest
				E		+= deltaE
				
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
	tDEFAULT			= 13   # Number of temperatures
	eDEFAULT			= 10   # Number of free-energy bins
	nDEFAULT			= 1000 # Number of sampled sequences
	parser = argparse.ArgumentParser(description='Samples random sequences from the uniform distribution and sorts them into -DDG/RT bins using a free-energy scoring matrix. The sequences are first sampled from the Boltzmann distribution e^{E/T} using Metropolis-Hastings sampling and then down-sampled to the uniform distribution using bin-specific rejection-sampling. Outputs the sampled sequences in FORMAT. The sequence identifiers contain energy-bin index iE, where 1 is the lowest-affinity bin, and a sequence index iSeq.')
	parser.add_argument('matrixFile', metavar='scoringMatrix.tsv', help='(Mononucleotide) scoring matrix. (COL 1) = "A,C,G,T", (COL 2..) = -ddG/RT values')
	parser.add_argument('-t', metavar='nTemp', help='Number of temperatures to use in the Metropolis-Hastings sampling (DEFAULT = %d)'%tDEFAULT, type=int, default=tDEFAULT)
	parser.add_argument('-e', metavar='nEnergyBins', help='Number of energy bins (DEFAULT = %d)'%eDEFAULT, type=int, default=eDEFAULT)
	parser.add_argument('-n', metavar='nSample', help='Number of sampled sequences (DEFAULT = %d)'%nDEFAULT, type=int, default=nDEFAULT)
#	parser.add_argument("--verbose", help="Increase output verbosity", action="store_true")
	args = parser.parse_args()

	#1. Reads and processes scoring matrix
	#1.1 - Reads matrix file 
	rawScoringMatrix	= sl.loadScoringMatrix(args.matrixFile)
	if rawScoringMatrix.keys() != [0]:
		sl.err("Model file should only contain scoring matrix (mononucleotides).")

	#1.2. Shifts and scales the scoring matrix so the sequence-score is in the interval [-1,0]
	scoringMatrix		= np.array([ i-max(i) for i in rawScoringMatrix[0] ]) / sum([ max(i)-min(i) for i in rawScoringMatrix[0] ])

	#2. Samples sequences at each temperature
	k            		= len(scoringMatrix) # Width of scoring matrix
	dE			= 1. / args.e        # Width of energy bins
	betaList 			= 2 * np.log(math.pow(4,k)) *  np.linspace(-1,1,args.t)      #List of inverse temperatures to loop over
	binnedSeq			= [ [ [] for ei in range(args.e) ] for ti in range(args.t) ] # (Empty) table of  energy-binned sequences

	for betaI in range(len(betaList)): #Loops over inverse temperatures.
		#2.1 Samples sequences with using Metropolis-Hastings sampling from distribution p[s] ~ e^(beta*E) 
		beta			= betaList[betaI]
		seqList			= sampleSequence( scoringMatrix, beta, 3*args.n )

		for s in seqList:
			#2.2 Identifies energy bin of sequence
			energyI 	= min(int(np.floor((s[1]+1)*args.e)), args.e-1)
			#2.3 Saves the sequence with probability proportional to e^(-beta*E). This downsamples to the uniform distribution within each bin.
			if (beta>0 and np.random.rand()<np.exp(-beta*( s[1]-(-1+energyI*dE)))) or (beta<=0 and np.random.rand()<np.exp(-beta*(s[1]-(-1+(energyI+1)*dE)))):
					binnedSeq[betaI][energyI] += [s]

	#3. Output
	#3.1 Converts from numeric to string representation
	outSeqs=[]
	for iE in range(args.e):
		outSeqs+=[[]]
		for Ti in range(args.t):
			for s in binnedSeq[Ti][iE]:
				outSeqs[iE]+=["".join(nucl[n] for n in s[0])]
		outSeqs[iE] =  np.random.permutation(outSeqs[iE])

	#3.2 Prints the sequences.
	for iE in range(args.e):
		for iS in range(min(args.n,len(outSeqs[iE]))):
			print ">iE=%d,iSeq=%d"%(iE+1,iS+1)
			print outSeqs[iE][iS]
			
##################  END #################

main()
