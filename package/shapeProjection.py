#!/usr/bin/env python

import numpy as np
import numpy.linalg as la
from scipy.optimize import minimize
import shapeLibrary as sl
import argparse

nucl = ("A","C","G","T")

##################  CLASS FOR EVALUATING VALUE AND GRADIENT OF MODEL PROJECTION  #################

class projectionModel:
	
	agnosticMono					= None 	# Agnostic model, free energy, mono-nucleotide part
	agnosticDi						= None	# Agnostic model, free energy, di-nucleotide part
	agnosticMonoAffinity			= None	# Scoring matrix exponentiated
	agnosticDiAffinity				= None	# Scoring matrix exponentiated
	diDesignMatrix					= None	# Design matrix for projected model
	L								= None	# Width of agnostic model
	k								= None	# Number of 
	useL1Shape 						= False # Activates L1 penalization of shape-sensitivity coefficients 
	useL1Mono  						= False	# Activates L1 penalization of mono
	lambdaL1Shape 					= 0 	# Lambda term for L1 penalization of shape
	lambdaL1Mono  					= 0		# Lambda term for L1 penalization of mono
	lambdaL2Mono					= 0		# Lambda term for L2 penalization of shape
	lambdaL2Shape					= 0		# Lambda term for L2 penalization of mono
	objectiveFunction				= None  # String, etiher 'KLDivergence', or 'affinity'
	gradMonoList1					= None	# Lists containing shape-constraint gradients, L1, Mono
	gradMonoList2					= None  # L2, Mono
	gradShapeList1					= None  # L1, Shape
	gradShapeList2					= None  # L2, Shape
	edgeReadout						= None	# Should readout of binding-site edge be included?
	iShape							= None	# Index of first shape sensitivity coefficient in parameter vector
	iMono							= None  # Same but mononucleotide
	iIntercept						= None  # Index of constant term
	iL1Shape						= None  # Index of the first dual variable for the L1 Shape penalty
	iL1Mono 						= None  # Index of the first dual variable for the L1 Mononucleotide penalty
	iEnd    						= None  # Last entry

	xMono = np.array([						# Array used to build projected mono matrix
		[ 3.,-1.,-1.,-1.],
		[-1., 3.,-1.,-1.],
		[-1.,-1., 3.,-1.]])/3.

	def __init__(self,agnosticModel, shapeModel, lambdaL1Mono=0., 
		lambdaL1Shape=0., lambdaL2Mono=0., lambdaL2Shape=0., objectiveFunction='KLDivergence',
		normalizeModel=True, edgeReadout=False):

		#Saves mechanism-agnostic (true) model
		if sorted(agnosticModel.keys()) != [0, 1]:	
			sl.err("Mechanism-agnostic model file must (only) have mono- and di-nucleotide values.")
		self.agnosticMono			= np.copy(agnosticModel[0])		#L*4 array
		self.agnosticMonoAffinity	= np.exp(self.agnosticMono)
		if len(agnosticModel[1].shape) == 3:
			self.agnosticDi			= np.copy(agnosticModel[1])
		else:
			self.agnosticDi			= np.array([ [ [ vi[4*n1+n2] for n2 in range(4) ] for n1 in range(4)] for vi in agnosticModel[1]])	#(L-1)*4*4 array
		
		self.agnosticDiAffinity		= np.exp(self.agnosticDi)
		self.L						= len(self.agnosticMono)
		self.edgeReadout			= edgeReadout


		#Sets the mean affinity to 1.0
		if normalizeModel:
			meanAffinity            =  self.computeSum(1, self.agnosticMonoAffinity, self.agnosticDiAffinity) / 4**self.L
			self.agnosticMono      -= np.log(meanAffinity) / self.L
			self.agnosticMonoAffinity = np.exp(self.agnosticMono)

		#Processes shape model
		if sorted(shapeModel.keys()) != [0, 1]:	
			sl.err("Shape model file must (only) have mono- and di-nucleotide values.")
		if len(shapeModel[0]) < 2:
			sl.err("Shape model must 2bp long or wider.")
		self.diDesignMatrix			= self.buildDiShapeDesignMatrix(shapeModel, self.L)
		self.k						= len(self.diDesignMatrix)

		self.objectiveFunction      = objectiveFunction

		# Setting up variables for L1 Penalization 
		if lambdaL1Mono == 0.:
			self.useL1Mono			= False
			self.lambaL1Mono		= None
		else:
			self.useL1Mono			= True	
			self.lambdaL1Mono		= lambdaL1Mono

		if lambdaL1Shape == 0.:
			self.useL1Shape			= False
			self.lambdaL1Shape		= None
		else:
			self.useL1Shape			= True	
			self.lambdaL1Shape		= lambdaL1Shape

		self.lambdaL2Mono 			= lambdaL2Mono
		self.lambdaL2Shape			= lambdaL2Shape

		#Setting up indices for accessing state vector
		self.iShape					= 0
		self.iMono					= self.iShape + self.k 
		self.iIntercept				= self.iMono + 3 * self.L
		self.iL1Shape				= self.iIntercept + 1
		if self.useL1Shape:	
			self.iL1Mono			= self.iL1Shape + self.k
		else:
			self.iL1Mono			= self.iL1Shape 
		if self.useL1Mono:
			self.iEnd				= self.iL1Mono + 4*self.L
		else:
			self.iEnd				= self.iL1Mono
	
		return

	def buildDiShapeDesignMatrix(self, shapeModel, L):
		""" Builds one dinucleotide scoring matrix (L*(L-1)*4*4) for each shape sensitivity coefficients (L)"""
		k             = len(shapeModel[0])

		#1. Writes the shape model using only dinucletide features
		shapeModelDi  = np.array( [ [ [ sml[4*n1+n2] for n2 in range(4) ] for n1 in range(4) ] for sml in shapeModel[1] ] )
		#1.1 Translates mononucleotide matrix (k*4) to dinucleotide
		monoDiRepTemp = [ np.array([sml for i in range(4)]) for sml in shapeModel[0] ]
		leftWeights   = [1.]+([0.5]*(k-2)) + [0.] 
		rightWeights  = [0.]+([0.5]*(k-2)) + [1.]
			
		shapeModelDi += np.array([ leftWeights[i]   *monoDiRepTemp[i].transpose() 
								+  rightWeights[i+1]*monoDiRepTemp[i+1]           for i in range(k-1)])

		
		#2. Creates one copy of the dinucleotide shape model for each position with a shape sensitivity coefficient 
		if self.edgeReadout:
			f = int((k-1)/2) # The shape model extend (and is cropped) (L-1)/2 basepairs beyond the mechanism-agnostic model
			diShapeDesignMatrix = np.array([ np.concatenate((np.zeros((x,4,4)), shapeModelDi, np.zeros((L+2*f-k-x,4,4))), axis=0)[f:-f] for x in range(L-k+2*f+1) ])
		else:
			f = 0
			diShapeDesignMatrix = np.array([ np.concatenate((np.zeros((x,4,4)), shapeModelDi, np.zeros((L+2*f-k-x,4,4))), axis=0)       for x in range(L-k+2*f+1) ])


		return diShapeDesignMatrix
	
	def vectorToMonoMatrix(self, v):
		""" Parses a state vector and returns a monucleotide scoring matrix """
		return np.array([ sum([self.xMono[i]*v[self.iMono+3*x+i] for i in range(3)]) for x in range(self.L) ])
	
	def vectorToDiMatrix(self, v):
		""" Parses a state vector and returns a dinucleotide scoring matrix """
		return sum([self.diDesignMatrix[i]*v[self.iShape+i] for i in range(self.k)]) 
	
	def vectorToModel(self, v):
		""" Parses a state vector and returns mono+di model"""
		return {0:self.vectorToMonoMatrix(v), 1:self.vectorToDiMatrix(v)}

	def computePartialSums(self, intercept, monoMatrix, diMatrix):
		""" Computes the left and right partial sums """
		(leftSum, rightSum)  = ([np.ones(4)* intercept], [np.ones(4)])
		
		for i in range(1, self.L):
			leftSum += [np.dot(leftSum[-1]*monoMatrix[i-1], diMatrix[i-1])]
			rightSum+= [np.dot(diMatrix[self.L-i-1], rightSum[-1]*monoMatrix[self.L-i])]

		return (np.array(leftSum), np.array(rightSum[::-1]))
	
	def computeSum(self, intercept, monoMatrix, diMatrix):
		""" Computes the sum over affinities"""
		partialSum = np.ones(4) * intercept
		for i in range(1, self.L):
			partialSum = np.dot(partialSum*monoMatrix[i-1], diMatrix[i-1])

		return np.dot(partialSum, monoMatrix[-1])

	def computeAffinityError(self, v, sign=1.0):
		""" Computes the affinity error (Not includig penalties) """

		projectedMono         = self.vectorToMonoMatrix(v)
		projectedDi           = self.vectorToDiMatrix(v)
		projectedMonoAffinity = np.exp(projectedMono)
	
		projectedDiAffinity   = np.exp(projectedDi)

		intercept11           = 1.
		mono11                = self.agnosticMonoAffinity**2 / 4 
		di11                  = self.agnosticDiAffinity**2

		intercept22           = np.exp(v[self.iIntercept]*2)  
		mono22                = projectedMonoAffinity**2      / 4
		di22				  = projectedDiAffinity**2

		intercept12           = np.exp(v[self.iIntercept])
		mono12				  = self.agnosticMonoAffinity * projectedMonoAffinity / 4
		di12                  = self.agnosticDiAffinity * projectedDiAffinity

		s11                   = self.computeSum(intercept11, mono11, di11)	
		s22                   = self.computeSum(intercept22, mono22, di22) 
		s12                   = self.computeSum(intercept12, mono12, di12)
		return (s11 + s22 - 2*s12) # / 4**self.L

	def computeAffinityErrorGradient(self, v, sign=1.0):
		""" Computes the gradient of the affinity error (Not includig penalties) """

		#Prepares afffinities
		projectedMono         = self.vectorToMonoMatrix(v)
		projectedDi           = self.vectorToDiMatrix(v)
		projectedMonoAffinity = np.exp(projectedMono)
		projectedDiAffinity   = np.exp(projectedDi)

		#Prepares affinity combinations
		intercept22           = np.exp(v[self.iIntercept]*2)
		mono22                = projectedMonoAffinity**2 / 4
		di22				  = projectedDiAffinity**2

		intercept12           = np.exp(v[self.iIntercept])
		mono12				  = self.agnosticMonoAffinity * projectedMonoAffinity / 4
		di12                  = self.agnosticDiAffinity * projectedDiAffinity

		#Computing partial sums
		(left22, right22)     = self.computePartialSums(intercept22, mono22, di22) 
		(left12, right12)     = self.computePartialSums(intercept12, mono12, di12)

		#Computes components intercept gradient
		rawInterceptGradient22 = np.dot(left22[-1], mono22[-1]) 
		rawInterceptGradient12 = np.dot(left12[-1], mono12[-1]) 

		#Computes components of mononucleotide gradient
		rawMonoGradient22     = left22 * right22 * mono22
		rawMonoGradient12     = left12 * right12 * mono12
	
		#Computes components of dinucleotide gradient
		rawDiGradient22       = np.array([ np.outer(left22[i]*mono22[i], right22[i+1]*mono22[i+1]) * di22[i] for i in range(self.L-1)])
		rawDiGradient12       = np.array([ np.outer(left12[i]*mono12[i], right12[i+1]*mono12[i+1]) * di12[i] for i in range(self.L-1)])

		#Combines components into full gradient
		rawInterceptGradient  = 2 * (rawInterceptGradient22 - rawInterceptGradient12 ) # / 4**self.L
		rawMonoGradient       = 2 * (rawMonoGradient22      - rawMonoGradient12      ) # / 4**self.L
		rawDiGradient         = 2 * (rawDiGradient22        - rawDiGradient12        ) # / 4**self.L

		#Computes Gradients in vector format
		vectorInterceptGradient = np.array([rawInterceptGradient])	
		vectorMonoGradient    = np.dot(rawMonoGradient, np.transpose(self.xMono)).flatten()
		vectorShapeGradient   = np.array([np.sum(rawDiGradient * self.diDesignMatrix[i]) for i in range(self.k)])

		out                   = np.zeros(1 + self.L*3 + self.k)
		out[self.iIntercept:self.iIntercept+1] = vectorInterceptGradient
		out[self.iMono:self.iMono+3*self.L]    = vectorMonoGradient
		out[self.iShape:self.iShape+self.k]    = vectorShapeGradient

		return out

	def computeKLDivergence(self, v, sign=1.0):
		""" Computes the KL Divergence (Not including penalties) """

		#Prepares afffinities
		projectedMono         = self.vectorToMonoMatrix(v)
		projectedDi           = self.vectorToDiMatrix(v)

		#Prepares affinity combinations
		intercept1            = 1.
		mono1				  = self.agnosticMonoAffinity / 4
		di1                   = self.agnosticDiAffinity 

		intercept2            = np.exp(v[self.iIntercept])
		mono2                 = np.exp(projectedMono) / 4
		di2				  	  = np.exp(projectedDi)
 
		#Computing partial sums
		(left1, right1)       = self.computePartialSums(intercept1, mono1, di1)
		(left2, right2)       = self.computePartialSums(intercept2, mono2, di2) 

		#Computes sums	
		s1                    = np.dot(left1[-1], mono1[-1])
		s2                    = np.dot(left2[-1], mono2[-1])

		#Computes components intercept gradient
		rawInterceptGradient1 = s1

		#Computes components of mononucleotide gradient
		rawMonoGradient1      = left1 * right1 * mono1
	
		#Computes components of dinucleotide gradient
		rawDiGradient1        = np.array([ np.outer(left1[i]*mono1[i], right1[i+1]*mono1[i+1]) * di1[i] for i in range(self.L-1)])

		#Combines components into full gradient
		out = 0.
		out +=         (0.                - v[self.iIntercept] ) * rawInterceptGradient1 / s1
		out += np.sum( (self.agnosticMono -  projectedMono     ) * rawMonoGradient1      / s1 )
		out += np.sum( (self.agnosticDi   -  projectedDi       ) * rawDiGradient1        / s1 )
		out += - np.log(s1/s2)
	
		return out

	def computeKLDivergenceGradient(self, v, sign=1.0):
		""" Computes the gradient of the KL Divergence (Not including penalties) """

		#Prepares afffinities
		projectedMono         = self.vectorToMonoMatrix(v)
		projectedDi           = self.vectorToDiMatrix(v)

		#Prepares affinity combinations
		intercept1            = 1.
		mono1				  = self.agnosticMonoAffinity / 4
		di1                   = self.agnosticDiAffinity 

		intercept2            = np.exp(v[self.iIntercept])
		mono2                 = np.exp(projectedMono) / 4
		di2				  	  = np.exp(projectedDi)
 
		#Computing partial sums
		(left1, right1)       = self.computePartialSums(intercept1, mono1, di1)
		(left2, right2)       = self.computePartialSums(intercept2, mono2, di2) 

		#Computes sums	
		s1                    = np.dot(left1[-1], mono1[-1])
		s2                    = np.dot(left2[-1], mono2[-1])

		#Computes components intercept gradient
		rawInterceptGradient1 = s1
		rawInterceptGradient2 = s2

		#Computes components of mononucleotide gradient
		rawMonoGradient1      = left1 * right1 * mono1
		rawMonoGradient2      = left2 * right2 * mono2
	
		#Computes components of dinucleotide gradient
		rawDiGradient1        = np.array([ np.outer(left1[i]*mono1[i], right1[i+1]*mono1[i+1]) * di1[i] for i in range(self.L-1)])
		rawDiGradient2        = np.array([ np.outer(left2[i]*mono2[i], right2[i+1]*mono2[i+1]) * di2[i] for i in range(self.L-1)])

		#Combines components into full gradient
	
		rawInterceptGradient  = rawInterceptGradient2 / s2 - rawInterceptGradient1 / s1 
		rawMonoGradient       = rawMonoGradient2      / s2 - rawMonoGradient1      / s1 
		rawDiGradient         = rawDiGradient2        / s2 - rawDiGradient1        / s1 

		#Computes Gradients in vector format
		vectorInterceptGradient = np.array([rawInterceptGradient])	
		vectorMonoGradient    = np.dot(rawMonoGradient, np.transpose(self.xMono)).flatten()
		vectorShapeGradient   = np.array([np.sum(rawDiGradient * self.diDesignMatrix[i]) for i in range(self.k)])

		out                   = np.zeros(1 + self.L*3 + self.k)
		out[self.iIntercept:self.iIntercept+1] = vectorInterceptGradient
		out[self.iMono:self.iMono+3*self.L]    = vectorMonoGradient
		out[self.iShape:self.iShape+self.k]    = vectorShapeGradient

		return out

	def func(self, v, sign=1.0):
		""" Computes objective function """

		out = 0.

		#Compute non-penalized objective fuction
		if self.objectiveFunction == 'affinityError':
			out += sign * self.computeAffinityError(v)
		elif self.objectiveFunction == 'KLDivergence':
			out += sign * self.computeKLDivergence(v)
		else:
			sl.err('Invalid objective function')

		#Adding L2 penalty:
		if self.lambdaL2Mono != 0:
			out += sign * self.lambdaL2Mono  * sum([ la.norm(mi,2)**2 for mi in self.vectorToMonoMatrix(v)])

		if self.lambdaL2Shape != 0:
			out += sign * self.lambdaL2Shape * la.norm(v[self.iShape:self.iShape+self.k],2)**2

		#Adding L1 penalty:
		if self.useL1Mono:
			out += sign * self.lambdaL1Mono  * np.sum(v[self.iL1Mono:self.iL1Mono+4*self.L])
		if self.useL1Shape:
			out += sign * self.lambdaL1Shape * np.sum(v[self.iL1Shape:self.iL1Shape+self.k])

		return out


	def func_deriv(self, v, sign=1.0):
		""" Derivative of objective function """

		gradient = np.zeros(self.iEnd)
		
		#Compute non-penalized objective fuction
		if self.objectiveFunction == 'affinityError':
			gradient[:self.k+3*self.L+1] = sign * self.computeAffinityErrorGradient(v)
		elif self.objectiveFunction == 'KLDivergence':
			gradient[:self.k+3*self.L+1] = sign * self.computeKLDivergenceGradient(v)
		else:
			sl.err('Invalid objective function')

		#Adding L2 contribution to gradient
		if self.lambdaL2Mono != 0:
			gradient[self.iMono:self.iMono+self.L*3] +=  2 * self.lambdaL2Mono * np.dot(self.vectorToMonoMatrix(v), self.xMono.transpose()).flatten()
		if self.lambdaL2Shape != 0:
			gradient[self.iShape:self.iShape+self.k] +=  2 * self.lambdaL2Shape * np.array(v[self.iShape:self.iShape+self.k])

		#Adding gradinets of L1 variables
		if self.useL1Mono:
			gradient[self.iL1Mono:self.iL1Mono+4*self.L] = np.ones(4*self.L)*self.lambdaL1Mono
		if self.useL1Shape:
			gradient[self.iL1Shape:self.iL1Shape+self.k] = np.ones(self.k)*self.lambdaL1Shape

		return gradient

	def getCons(self, b=10):
		""" Generates constraints used to implement dual of L1 penalty """

		cons = []
		if self.useL1Mono:
			self.gradMonoList1 = []
			self.gradMonoList2 = []
			for i in range(self.L):

				for n in range(4):

					# L1Mono -u <= x  <=> u+x >= 0
					grad                                    = np.array([0.] * self.iEnd)
					grad[self.iMono+3*i:self.iMono+3*(i+1)] = self.xMono[:,n]
					grad[self.iL1Mono+4*i+n]                = 1.
					self.gradMonoList1+=[grad[:]]

					cons += [{'type': 'ineq', # -u<=x   <=> u+x >=0
    		        	 'fun' : lambda x,i=i,n=n: np.array([  np.dot(x[self.iMono+3*i:self.iMono+3*(i+1)], self.xMono[:,n] )  + x[self.iL1Mono+4*i+n] ]),
	        		     'jac' : lambda x,i=i,n=n: self.gradMonoList1[4*i+n]}]	

					# L1Mono x <= u   <=> u-x >= 0
					grad                                    = np.array([0.] * self.iEnd)
					grad[self.iMono+3*i:self.iMono+3*(i+1)] = -self.xMono[:,n]
					grad[self.iL1Mono+4*i+n]                = 1.
					self.gradMonoList2+=[grad[:]]

					cons += [{'type': 'ineq', # -u<=x   <=> u+x >=0
    		         	 'fun' : lambda x,i=i,n=n: np.array([ -np.dot(x[self.iMono+3*i:self.iMono+3*(i+1)], self.xMono[:,n] )  + x[self.iL1Mono+4*i+n] ]),
	        		     'jac' : lambda x,i=i,n=n: self.gradMonoList2[4*i+n]}]	

		if self.useL1Shape:
			self.gradShapeList1 = []
			self.gradShapeList2 = []
			for i in range(self.k):

				# L1Shape -u <= x  <=> u+x >= 0
				grad                  = np.array([0.] * self.iEnd)
				grad[self.iShape+i]   = 1.
				grad[self.iL1Shape+i] = 1.
				self.gradShapeList1+=[grad[:]]

				cons += [{'type': 'ineq', # -u<=x   <=> u+x >=0
    		         'fun' : lambda x,i=i: np.array([  x[self.iShape+i] + x[self.iL1Shape+i] ]),
        		     'jac' : lambda x,i=i: self.gradShapeList1[i]}]	

				# L1Shape x <= u   <=> u-x >= 0
				grad                  = np.array([0.] * self.iEnd)
				grad[self.iShape+i]   = -1.
				grad[self.iL1Shape+i] = 1.
				self.gradShapeList2+=[grad[:]]
				cons += [{'type': 'ineq', # -u<=x   <=> u+x >=0
    		         'fun' : lambda x,i=i: np.array([ -x[self.iShape+i] + x[self.iL1Shape+i] ]),
        		     'jac' : lambda x,i=i: self.gradShapeList2[i]}]	

		return tuple(cons)

	def scoreSeq(self, inSeq, mutList=()):
		"""Scores a sequence using the mechanism-agnostic model"""
		seq = inSeq[:]
		for i in range(len(mutList)):
			seq[mutList[i][0]] = mutList[i][1]

		return sum([self.agnosticMono[i,seq[i]] for i in range(len(seq))]) + sum([self.agnosticDi[i,seq[i],seq[i+1]] for i in range(len(seq)-1)])
		
	def getMaxSeq(self, monoMatrix, diMatrix):
		""" Identifies the max-affinity sequence using dynamic programming"""
		state = [ [monoMatrix[0,i], [i]] for i in range(4) ]
		for x in range(0,self.L-1):
			state = [sorted([ [state[i][0]+diMatrix[x,i,j] +monoMatrix[x+1,j], state[i][1]+[j]  ] for i in range(4)], key=lambda x: x[0])[-1] for j in range(4)]
		maxSeq = sorted(state,  key=lambda x: x[0])[-1]
		return maxSeq

	def getSeed(self):
		""" Generates a naive seed"""

		#1) Identify top-scoring sequence
		maxSeq = self.getMaxSeq(self.agnosticMono, self.agnosticDi)
		
		#2) Compute matrix of single-mutation deviations.
		mutMatrix =  [[self.scoreSeq(maxSeq[1], ((x,i),) )-maxSeq[0] for i in range(4)] for x in range(self.L)]
		
		#3) Phi1 = 0.75 * {E[A]-E[T], E[C]-E[T], E[G]-E[T]}
		monoParameters =  [ 0.75 * ( mutMatrix[x][n]-mutMatrix[x][3] )  for x in range(self.L) for n in range(3)]
	
		#4) Phi0 = MaxSeq + Sum[T]
		interceptValue =  maxSeq[0] + sum( [ 0.25*sum(mutMatrix[x]) for x in range(self.L)])

		#Constructs output Vector
		outVector = [0.] * self.iEnd
		outVector[self.iMono:self.iMono+3*self.L] = monoParameters
		outVector[self.iIntercept] = interceptValue
		if self.useL1Mono:
			outVector[self.iL1Mono:self.iL1Mono+self.L] = [ la.norm(np.array(mi) - np.mean(np.array(mi)),1) for mi in mutMatrix ]

		return outVector
			
#	def getZeroSeed(self):
#		""" Generates seed to start the optimization at the origin"""
#		return [0. for i in range(self.iEnd)]

	#Naive bounds
	def getBounds(self, b=10):
		""" Generates bounds for the optimization """ 
		out = []
		for i in range(self.iEnd):
			if i == self.iIntercept:
				out += [(-3*b,3*b)]
			else:
				out += [(-b,b)]

		return tuple(out)

	def fitModel(self, n=1000, b=10, verbose=False, useCOBYLA=False):
		""" Performs the optimization in the shape projection."""
		if useCOBYLA:
			sl.disp("> Fitting using COBYLA.", verbose)
			res = minimize(self.func, self.getSeed(), args=(1.0,),
	               constraints=self.getCons(b=b), method='COBYLA', options={'disp': verbose, 'maxiter': n})
		else:
			sl.disp("> Fitting using SLSQP.", verbose)
			res = minimize(self.func, self.getSeed(), args=(1.0,), jac=self.func_deriv,
   		           constraints=self.getCons(b=b), method='SLSQP', options={'disp': verbose, 'maxiter': n}, bounds=self.getBounds(b)) 

		return res

	def formatModelCSV(self, v, reportL1=False):
		""" Formats the projected model as a string"""
		out  = []
		#Adding shape
		out += list(v[self.iShape:self.iShape+self.k])
		#Adding mono
		for i in range(self.L):
			out += list(np.dot(v[self.iMono+3*i:self.iMono+3*(i+1)], self.xMono))
		#Adding intercept 
		out += [v[self.iIntercept]]
		#Adding L1 norm variables
		if reportL1:
			if self.useL1Shape:
				out += list(v[self.iL1Shape:self.iL1Shape+self.k])
			if self.useL1Mono:
				out += list(v[self.iL1Mono:self.iL1Mono+self.L*4])
		
		return ",".join([ "%f"%v for v in out])
		
	def formatHeaderCSV(self, reportL1=False):
		""" Formats a header for the output"""
		out  = []
		#Adding shape
		out += [ "shape:%d"%(i+1) for i in range(self.k)]
		#Adding mono
		out += [ n+":%d"%(i+1) for i in range(self.L) for n in nucl]
		#Adding intercept
		out += ["intercept"]
		#Adding L1 norm variables
		if reportL1:
			if self.useL1Shape:
				out += [ "eta_shape:%d"%i for i in range(self.k)]
			if self.useL1Mono:
				out += [ "eta_"+n+":%d"%(i+1) for i in range(self.L) for n in nucl]
		
		return ",".join(out)


	def formatModelTable(self, v, reportL1=False):
		monoMatrix			= self.vectorToMonoMatrix(v).transpose()
		outStr 				= ""
		for ni in range(4):
			outStr			+= "%s\t%s\n"%(nucl[ni], "\t".join([ "%f"%ri for ri in monoMatrix[ni]]))
		outStr				+= "shape\t%s\n"%"\t".join([ "%f"%xi for xi in v[self.iShape:self.iMono]])
		if reportL1:
			if self.useL1Mono:
				for n in range(4):
					outStr	+= "eta_%s\t%s\n"%(nucl[n],"\t".join([ "%f"%v[self.iL1Mono+4*i+n] for i in range(self.L)]))
			if self.useL1Shape:
				outStr 		+= "eta_shape\t%s\n"%("\t".join([ "%f"%v[self.iL1Shape+i] for i in range(len(self.diDesignMatrix))]))

		return outStr.rstrip('\n')
   
##################  VARIOUS FUNCTIONS  #################

##################  MAIN FUNCTION  #################

def main():

	#DEFAULT values of parameters 
	nDEFAULT=1000
	bDEFAULT=10.
	(l1sDEFAULT, l1mDEFAULT, l2sDEFAULT, l2mDEFAULT)=("0.0", "0.0", "0.0", "0.0")
	#Creating parser
	parser = argparse.ArgumentParser(description='Performs shape projection. Uses L2 penalty terms by default.')
	parser.add_argument('agnosticFile', metavar='bindingModel.tsv', help='File containing mechanism-agnistic  mono- and di-nucleotide binding model.')
	parser.add_argument('shapeFile', metavar='shapeModel.tsv', help="File containing mono- and di-nucleotide model of shape paramter. Multiple files can be separated by ','")
	parser.add_argument('-n', metavar='nIt', help='Maximum number of iterations(DEFAULT = %d)'%nDEFAULT, type=int, default=nDEFAULT)
	parser.add_argument('-b', metavar='bound', help='Bound on parameter values (DEFAULT = %f)'%bDEFAULT, type=float, default=bDEFAULT)

	parser.add_argument('-L1Shape', metavar='lambda', help='Lambda parameter for shape readout L1 penalty. (DEFAULT = %s)'%l1sDEFAULT, type=str, default=l1sDEFAULT)
	parser.add_argument('-L1Mono', metavar='lambda', help='Lambda parameter for mono readout L1 penalty. (DEFAULT = %s)'%l1mDEFAULT,   type=str, default=l1mDEFAULT)
	parser.add_argument('-L2Shape', metavar='lambda', help='Lambda parameter for shape readout L2 penalty. (DEFAULT = %s)'%l2sDEFAULT, type=str, default=l2sDEFAULT)
	parser.add_argument('-L2Mono', metavar='lambda', help='Lambda parameter for mono readout L2 penalty. (DEFAULT = %s)'%l2mDEFAULT,   type=str, default=l2mDEFAULT)

	parser.add_argument("--edgeReadout", help="Includes shape-readout of first and last bps of binding site by truncating the sequence model.", action="store_true")
	parser.add_argument("--affinityError", help="Minimizes the KL Divergence (instead of the KL divergence)", action="store_true")
	parser.add_argument("--COBYLA", help="Uses the COBYLA minimizer algorithm (SLSQP default)", action="store_true")
	parser.add_argument("--header", help="Prints a header", action="store_true")
	parser.add_argument("--reportL1", help="Report the dual L1 variables", action="store_true")
	parser.add_argument("--verbose", help="Increase output verbosity", action="store_true")
	args = parser.parse_args()

	#Building list of penalties to sweep over
	if args.L1Shape=="0.0" and args.L2Shape=="0.0" and args.L2Shape=="0.0" and args.L2Mono=="0.0":
		lambdaL1ShapeList = [0.]
		lambdaL1MonoList  = [0.]
		lambdaL2ShapeList = [1.]
		lambdaL2MonoList  = [1.]
	else:
		lambdaL1ShapeList = [ float(d) for d in args.L1Shape.split(",")] 
		lambdaL1MonoList  = [ float(d) for d in args.L1Mono.split(",")] 
		lambdaL2ShapeList = [ float(d) for d in args.L2Shape.split(",")] 
		lambdaL2MonoList  = [ float(d) for d in args.L2Mono.split(",")] 


	#Building model
	agnosticModel	  = sl.loadScoringMatrix(args.agnosticFile)
	shapeModel        = sl.loadScoringMatrix(args.shapeFile)

	if args.affinityError:
		of = 'affinityError'
	else:
		of = 'KLDivergence'

	sl.disp("> Running initial fit to get scale of parameters.", args.verbose)
	model				= projectionModel(agnosticModel, shapeModel, objectiveFunction=of, edgeReadout=args.edgeReadout)
	res					= model.fitModel(n=args.n, b=args.b, verbose=args.verbose, useCOBYLA=args.COBYLA)
  
	#Determines the scaling factor for the L1 and L2 penalties
	unpenalizedValue	= model.func(res.x) 															#Function value of unpenalized model
	unpenalizedL2Mono 	= sum([ la.norm(mi,2)**2 for mi in model.vectorToMonoMatrix(res.x)])			#L2 norm of mono in unpenalized model
	unpenalizedL2Shape	= la.norm(res.x[model.iShape:model.iShape+model.k],2)**2						#L2 norm of shape in unpenalized model
	unpenalizedL1Mono 	= sum([ la.norm(mi,1)    for mi in model.vectorToMonoMatrix(res.x)])			#L1 norm of mono in unpenalized model
	unpenalizedL1Shape	= la.norm(res.x[model.iShape:model.iShape+model.k],1)							#L1 norm of shape in unpenalized model

	lambdaL1MonoList 	= list(4 * np.array(lambdaL1MonoList)  * unpenalizedValue / unpenalizedL1Mono)  #Lists of lambda values to try
	lambdaL2MonoList 	= list(4 * np.array(lambdaL2MonoList)  * unpenalizedValue / unpenalizedL2Mono)
	lambdaL1ShapeList 	= list(    np.array(lambdaL1ShapeList) * unpenalizedValue / unpenalizedL1Shape)
	lambdaL2ShapeList 	= list(    np.array(lambdaL2ShapeList) * unpenalizedValue / unpenalizedL2Shape)

	firstFit=True
	sl.disp("> Starting to fit.", args.verbose)

	#Prints in CSV format if the program runs and outputs multiple fits.
	printCSV			= len(lambdaL1ShapeList)*len(lambdaL1MonoList)*len(lambdaL2ShapeList)*len(lambdaL2MonoList) > 1

	for lambdaL1Shape in lambdaL1ShapeList: #Loops over lambda values	
		for lambdaL1Mono  in lambdaL1MonoList:	
			for lambdaL2Shape in lambdaL2ShapeList:	
				for lambdaL2Mono  in lambdaL2MonoList:	

					sl.disp(">> LAMBDA: (L1Shape, L1Mono, L2Shape, L2Mono) = (%f, %f, %f, %f)"%(lambdaL1Shape, lambdaL1Mono, lambdaL2Shape, lambdaL2Mono), args.verbose)

					model				= projectionModel(agnosticModel, shapeModel, lambdaL1Mono=lambdaL1Mono, lambdaL1Shape=lambdaL1Shape, lambdaL2Mono=lambdaL2Mono, lambdaL2Shape=lambdaL2Shape, objectiveFunction=of, edgeReadout=args.edgeReadout)

					if args.header and firstFit and printCSV:
						print "lambdaL1Shape,lambdaL1Mono,lambdaL2Shape,lambdaL2Mono,affinityError,KLDivergence,"+model.formatHeaderCSV(reportL1=args.reportL1)
						firstFit		= False

					res					= model.fitModel(n=args.n, b=args.b, verbose=args.verbose, useCOBYLA=args.COBYLA)

					outVector			= res.x

					#The intercept is unconstrained when KL divergence is minimized. Sets it using naive seed on projected model 
					if not args.affinityError:
						projectedModel  = projectionModel(model.vectorToModel(outVector), shapeModel, edgeReadout=args.edgeReadout)
						outVector[model.iIntercept] = projectedModel.getSeed()[projectedModel.iIntercept]

					affinityError		= model.computeAffinityError(outVector)
					klDivergence		= model.computeKLDivergence(outVector)

					if args.verbose:
						print res

					if printCSV:
						print "%f,%f,%f,%f,%f,%f,"%(lambdaL1Shape, lambdaL1Mono, lambdaL2Shape, lambdaL2Mono, affinityError, klDivergence)+ model.formatModelCSV(outVector, reportL1=args.reportL1)
					else:
						print model.formatModelTable(res.x, reportL1=args.reportL1)
				
	
##################  END #################

main()
