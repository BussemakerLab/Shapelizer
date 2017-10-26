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
def kMerToVector(kMer, mono, di, all, centralMono):
	"""Converts a kMer to a line in a design matrix"""

	
	if mono or all or centralMono:

		
		#Computes mononucleotide representation
		monoRep			= dict([ (nucl[i], np.identity(4)[i].tolist()) for i in range(len(nucl))] )
		x				= []

		#Determines what indices to use
		if centralMono:
			k           = len(kMer)
			iStart      = (k-1)/2
			iRange      = range(iStart,k-iStart)
		else:
			iRange      = range(len(kMer))
			
		#Builds x-vector
		for i in iRange:
			x			+=  monoRep[kMer[i]]

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
	parser = argparse.ArgumentParser(description='Uses linear regression to analyze a pentamer table.')
	parser.add_argument('kmerFile', metavar='kmerValue.csv', help='kMer file. (COL 1) = kmer, (COL 2) = values')
	parser.add_argument("--header", help="First line in kmer file is header.", action="store_true")
	parser.add_argument("--verbose", help="Increase output verbosity", action="store_true")
	paramGroup = parser.add_mutually_exclusive_group(required=True)
	paramGroup.add_argument('--mono', help="Use mononucleotides as predictors.", action='store_true')
	paramGroup.add_argument('--di', help="Use dinucleotides as predictors.", action='store_true')
	paramGroup.add_argument('--all', help="Use all base-base pairs as predictors.", action='store_true')
	paramGroup.add_argument('--centralMono', help="Use only central mononucleotides as predictors.", action='store_true')
	outGroup = parser.add_mutually_exclusive_group(required=True)
	outGroup.add_argument('--R2', help="Reports R^2", action='store_true')
	outGroup.add_argument('--crossR2', help="Reports the hold-one-out R^2", action='store_true')
	outGroup.add_argument('--betas', help="Reports the regression coefficients beta=(X^T*X)^-1*(X^T*y).", action='store_true')
	outGroup.add_argument('--yHat', help="Reports the predicted values X^T*beta", action='store_true')
	outGroup.add_argument('--residual', help="Reports the residual", action='store_true')
	args = parser.parse_args()

	(kMers, values, k, nCol, header) = sl.readKMerTable(args.kmerFile, args.header)

	if not args.crossR2: #(Not used when cross-validated R2 is computed...)

		#Performing initial mononucleotide-only regression
		if args.centralMono:
			xMono			= np.array([ kMerToVector(km, False, False, False, True ) for km in kMers] )
		else:
			xMono			= np.array([ kMerToVector(km, True,  False, False, False) for km in kMers] )

		#Performs linear regression: beta = Pseudoinverse(X^T*X)*X^T*X
		betaMono 			= sp.linalg.pinvh(xMono.transpose().dot(xMono)).dot(xMono.transpose().dot(values))
		#Computes the predicted value yHat = X^T*beta
		yHatMono			= xMono.dot(betaMono)

		#Performs secondary regression using interactions
		if args.di or args.all:
			xInt			= np.array([ kMerToVector(km, False, args.di, args.all, False) for km in kMers] )
			#Computing the number of parameters 
			nParam			= len([ np.abs(ev) for ev in np.linalg.eigvals(xInt.transpose().dot(xInt)) if np.abs(ev) > 1e-8 ])
			#Performs linear regression: beta = Pseudoinverse(X^T*X)*X^T*X
			betaInt 		= sp.linalg.pinvh(xInt.transpose().dot(xInt)).dot(xInt.transpose().dot(values-yHatMono))
			#Computes the predicted value yHat = X^T*beta
			yHatInt			= xInt.dot(betaInt)
			#Full predicted value
			yHatFull		= yHatMono + yHatInt
		else:
			nParam			= len([ np.abs(ev) for ev in np.linalg.eigvals(xMono.transpose().dot(xMono)) if np.abs(ev) > 1e-8 ])
			yHatFull		= yHatMono

	# WRITES OUTPUT
	###############
	if args.R2:
		# PRINTS CORRELATION COEFFICIENTS 
		#################################
		#Prints the R^2
		if args.header:
			print ",".join(header[1:])
		print ",".join(["%f"%math.pow(st.pearsonr(yHatFull[:,i], values[:,i])[0],2) for i in range(len(values[0]))])

	elif args.crossR2:
		# COMPUTES AND PRINTS THE HOLD-ONE-OUT CROSS-VALIDATED R2
		#########################################################
		xAll				=  np.array([ kMerToVector(km, args.mono, args.di, args.all, args.centralMono) for km in kMers] )

		#Identifies reverse-complement pairs of sequences that shouldbe held out 
		kmerPairs			= {}
		trantab				= maketrans("ACGT", "TGCA")
		for i in range(len(kMers)):
			km				= kMers[i]
			rcKm			= km[::-1].translate(trantab)
			if rcKm not in kmerPairs:
				kmerPairs[kMers[i]]	= (i, kMers.index(rcKm))

		#Loops over kmers pair and holds out one at the time
		inValues			= []
		outValues			= []
		for pair in kmerPairs.values():
			xTemp 			= np.delete(xAll,   pair, 0)
			valuesTemp		= np.delete(values, pair, 0)

			inValues		+= values[pair,].tolist()
			outValues 		+= (xAll[pair,:].dot(sp.linalg.pinvh(xTemp.transpose().dot(xTemp)).dot(xTemp.transpose().dot(valuesTemp)))).tolist()

		#Prints mean r^2:
		if args.header:
			print ",".join(header[1:])
		print ",".join(["%f"%pow(st.pearsonr(np.array(inValues)[:,i], np.array(outValues)[:,i])[0],2) for i in range(nCol) ])

	elif args.betas:
		#PRINTS REGRESSION COEFFICIENTS
		###############################
		if len(values[0])>1:
			sl.err("--betas option only possible when k-mer table has one data column")
		
		diNucl				= [ n1+n2 for n1 in nucl for n2 in nucl ]

		#Organizes coefficients.
		model				= {0:{}}
		#Adds data from mononucleotide coefficients
		for i in range(len(nucl)):
			model[0][nucl[i]] 		= [ betaMono[4*j+i,0] for j in range(k)]

		if args.di:
			#Adds dinucleotide entries
			model[1] 				= {}
			for i in range(len(diNucl)):
				model[1][diNucl[i]] = [ betaInt[16*j+i,0] for j in range(k-1)]

		elif args.all:
			#Adds interacton entries
			for d in range(1,k):
				model[d] 			= {}
				for di in diNucl:
					model[d][di]	= [0.] * (k-d)
			
			#Fills in values in model
			for i1 in range(k):
				for i2 in range(k):
					for n1 in range(len(nucl)):
						for n2 in range(len(nucl)):
							bi										= betaInt[(4*k)*(4*i1 + n1) + 4*i2 + n2, 0]
							if i1 == i2 and n1 == n2:
								#This should be zero if all k-mers are present in the k-mer table
								model[0][nucl[n1]][i]				+= bi
							elif i2>i1:
								model[i2-i1][nucl[n1]+nucl[n2]][i1] += bi
							elif i2<i1:
								model[i1-i2][nucl[n2]+nucl[n1]][i2] += bi

		#Moves mean value to the intercept
		m 							= np.mean(values[:,0])
		model["intercept"] 			= m
		for n in nucl:
			for i in range(k):
				#Moves mean to intercept
				model[0][n][i] 		-= m/k

		#Prints model
		print "intercept\t%f"%model["intercept"]
		#Mono
		for n in nucl:
			print n+"\t"+("\t".join(["%f"%v for v in model[0][n]]))
		#Di
		for d in range(1,k):
			if d in model:
				for di in diNucl:
					print di+":%d\t"%d+("\t".join(["%f"%v for v in model[d][di]]))
				
	elif args.yHat:

		# PRINTS THE PREDICTED VALUES
		#############################
		if args.header:
			print ",".join(header)
		for i in range(len(kMers)):
			print kMers[i]+","+(",".join(["%f"%yi for yi in yHatFull[i]]))

	elif args.residual:

		# PRINTS THE RESIDUAL VALUES
		############################
		if args.header:
			print ",".join(header)
		for i in range(len(kMers)):
			print kMers[i]+","+(",".join(["%f"%(values[i,j]-yHatFull[i,j]) for j in range(len(yHatFull[i]))]))
		
##################  END #################

main()
