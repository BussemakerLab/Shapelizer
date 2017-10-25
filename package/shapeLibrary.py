import sys
import numpy as np

nucl = ("A", "C", "G", "T")

#Ouputs an error message to STDERR and quits
def err(msg):
	sys.stderr.write("ERROR: %s\n"%msg)
	exit(1)

#Ouputs a message to STDERR if verbos argument is set to true(Default)
def disp(msg, verbose=True):
	if verbose:
		sys.stderr.write(msg+"\n")

#UPDATE TO RESPECT THE CURRENT FORMAT
def loadScoringMatrix(path, nucl=["A", "C", "G", "T"], readIntercept=False):
	""" Reads a scoring matrix file. Has format:
A	1	2	-1
C	0	2	1
G	-1	-1	2
T	-2	-2	4
AA:1	1	3
..."""

	nucl   = ["A", "C", "G", "T"]	
	dinucl = [ n1+n2 for n1 in nucl for n2 in nucl]

	#Reads data file, organizes values by 1) interaction distance, 2) bases.
	organizedData = {}
	intercept     = None
	with open(path) as f:
		for line in f:
			columnValues  = line.rstrip().split('\t')
			splitKey       = columnValues[0].split(":")
			if len(splitKey) == 1:
				if columnValues[0].lower() == 'intercept':
					if readIntercept:
						intercept = float(columnValues[1])
					continue
				n = columnValues[0]
				d = 0
			else:
				n = splitKey[0]
				d = int(splitKey[1])
			
			if d not in organizedData:
				organizedData[d] = {}

			organizedData[d][n] = [ float(x) for x in columnValues[1:]]

	#Organized the values into matrices
	outData = {}
	for d in organizedData:
		if d==0:
			if sorted(organizedData[d].keys()) != nucl:
				print (organizedData[d].keys())	
				print nucl
				err("Scoring matrix file contains incomplete set of mononucleotides")
			outData[d] = np.array([ organizedData[d][n] for n in nucl ]).transpose()
		else:
			if sorted(organizedData[d].keys()) != dinucl:
				err("Scoring matrix file contains incomplete set of dinucleotides")
			outData[d] = np.array([ organizedData[d][n] for n in dinucl ]).transpose()
			
	if readIntercept:
		if intercept is None:
			outData['intercept'] = 0.0
		else:
			outData['intercept'] = intercept
	return outData




def readKMerTable(filePath, hasHeader=False):
	"""Parses kmer file"""
	header				= None
	k					= None
	nCol				= None
	kMers				= []
	values 				= []

	if filePath == "-":
		f = sys.stdin
	else:
		f = open(filePath)


	if hasHeader:
		header		= f.readline().rstrip().split(',')
			
	for l in f:
		d			= l.rstrip().split(",")
		km			= d[0].strip('"').upper()

		if k is None:
			k		= len(km)
			nCol	= len(d)-1
		elif k != len(km):
			err('All kmers must be of equal length')

		kMers		+= [km]
		try:
			values  += [np.array([ float(di) for di in d[1:] ])]
		except ValueError:
			err("Table entries must be float numbers. LINE='%s'."%(l.rstrip()))

	if f is not "-":
		f.close()

	values				= np.array(values)

	return (kMers, values, k, nCol, header)

def printMatrix(matrix, f=sys.stdout):
	for i in range(len(matrix)):
		if i==0:
			f.write("[[")
		else:
			f.write(" [")
		f.write(", ".join([ "%f"%cv for cv in matrix[i] ]))
		if i==len(matrix)-1:
			f.write("]]\n")
		else:
			f.write("],\n")

def scoreSeq(seq, model):
	out = 0
	#Adds intercept
	if 'intercept' in model:
		out += model['intercept']

	#Scores mononucleotide matrix
	baseRep = dict([ (nucl[i], i) for i in range(4) ])
	k     =  len(model[0])
	for x in range(k):
		out += model[0][x, baseRep[seq[x]]]
	
	#Scores base-base interactions
	dList = [ di for di in model.keys() if di != 'intercept' and di > 0]
	for d in dList:
		for x in range(k-d):
			out += model[d][x][4*baseRep[seq[x]]+baseRep[seq[x+d]]]

	return out

