from Bio import SeqIO, AlignIO
from generateTree import loadTree, loadSequences, getSequencesOf, getLineagesOf
from helpers import *
from itertools import chain, izip, combinations

import logging
import numpy as np
import os
import random
import re
import subprocess
import sys



def replaceInnerGap(line):
	'''Replaces the internal gaps '-' of a string by ':'.

	 	line:	"" A sequence to edit.
	 
	 	return	"" The string with inner gaps replaced by ':'.
	'''
	frstA = re.search(r'^-*[^-]', line).end() - 1
	lastA = re.search(r'[^-]-*$', line).start()
	return line[ :frstA] + (line[frstA:lastA].replace('-', ':')) \
		+ line[lastA: ]


def computeDist_outterGapConserved(aligns):
	'''Computes the distance between the sequences in a .aln file,
		counting padding gaps as wildcard matches.

	 	aligns	 	[""] List of seqs to align.
	 
	 	returns 	The mismatch % between aligned sequences.
	'''
	miss = 0
	succ = 0
	for i in range(len(aligns)):
		aligns[i] = replaceInnerGap(aligns[i]) 

	for col in izip(*aligns):
		colSet = set(col)
		setSize = len(colSet)
		if ':' in colSet or (setSize > 2) or (setSize == 2 and not '-'\
			 in colSet):
			miss = miss + 1
		else:
			succ = succ + 1
	return miss / (miss + float(succ))


def computeDist(aligns, distFn):
	'''Computes the distance between sequences in a .aln file, using the 
		distFn specified.

	 	aligns: 	[] List of seqs to align""
	 	distFn: 	"" Uppercase ID for the dist function to use
	 	
	 	return		The distance as a % between sequences in an 
					alignFile
	'''
	distFns = {"OUTTER_GAP_CONSERVED": computeDist_outterGapConserved}
	
	if distFns.has_key(distFn):
		return distFns[distFn](aligns)
	else:
		raise NameError("Error: Invalid grading Fn.")


def computePairwiseDistStats(alignedSeqs, distFn):
	'''Does pairwise subsampling on a list of sequences.  Prints the 
		distnace matrix and computes summary stats. 

	 	alignedSeqs:	[] The sequences to compare.

		return		(,) A tuple of (min,max,avg,std) distances
					summarizing the table.
	'''
	sampleSize = len(alignedSeqs)
	distMatrix = np.zeros((sampleSize, sampleSize))
	valList = []
	for i,j in combinations(range(sampleSize), 2):
		distMatrix[i,j] = computeDist([alignedSeqs[i],alignedSeqs[j]],
						distFn)
		valList.append(distMatrix[i,j])

	distMatrix
	myMin = min(valList)
	myMax = max(valList)
	myAvg = np.mean(valList)
	myStd = np.std(valList)
	sol = {"min": myMin, "max":myMax, "avg":myAvg, "std":myStd}
	
	logging.debug("Pairwise Data:\n")
	for key in sol.keys():
		logging.debug("%s: %f\n" % (key, sol[key]))
	return sol


def emptySampleError(varName):
	logging.error("%s must be greater than or equal to 0." % varName)
	sys.exit()

def computeAlnDistance(t, s, tax, upToLvl, refDB, distFn, outDir="tmp", 
				maxPoolSize=100, pairwiseSampleSize=25):
	'''Computes the min, max, avg % genetic variance amongst members of a \
		given sequences.  Primary function.

	 	t		{{}{}} A tree from loadTree()
	 	s		{{}{}} A tree form loadSequences()
	 	refDB 		{} A SeqIO fast index.  Result of calling 
					indexFasta()
	 	tax		"a;b;c" A semicolon-delimited taxonomy string
	 	upToLvl		# A digit between 1 and 8 indicating depth of 
					tax
	 	distFn		"" Uppercase ID for the dist function to use
		outDir		"" Name of the local directory to use for file
					 output.	 	
		maxPoolSize	# A number indicating the maxmimum size of a
					 the taxa to multi-seq align.
	 	pairwiseSampleSize	# A number indicating the maxmimum 
						size of a subsample from the 
						tax pool to pariwise-seq align.

	 	returns		#.## Percentage diversity between sequences 
					at a given taxanomic level
	'''
	# poolSize must be positive or this is all pointless
	if maxPoolSize < 2:
		emptySampleError("maxPoolSize")
		os.exit()
	# same goes for pairwiseSampleSize
	if pairwiseSampleSize < 2:
		emptySampleError("subSampleSize")
		os.exit()

	# A safe file name delimiter to use as a replacement for unsafe delims.
	fileNameDelim = "_"
	# The directory to write to
	ensureDirExists(outDir)
	hitsFastaFilePath = "%s/%s.fasta" % (outDir, sanitizeFileName(tax, "_"))
	
	child_seqs = chain.from_iterable(getSequencesOf(t, s, tax, upToLvl))
	logging.info("Fetching sequences from database...")
	descendants = list(searchFastaByID(refDB, child_seqs))
	numDescendants = len(descendants)
	if numDescendants < 2:
		logging.warning("Ignoring a lineage of size %d: %s." % 
							(numDescendants, tax))
		return(tax,-1,-1,-1,-1)
	
	if numDescendants < maxPoolSize:
		logging.warning("The number of available taxonomic\
			matches is smaller than maxPoolSize.  \
			Ignoring maxPoolSize.")

	samplePool = randomSubSample(descendants, maxPoolSize)

	# publish matching seqs
	logging.info("Writing results to %s ..." % hitsFastaFilePath)
	SeqIO.write(samplePool, open(hitsFastaFilePath, 'w'), "fasta")

	# align seqs
	alignFile = muscleAlign(hitsFastaFilePath)
	aligns = [x for x in parseAln(alignFile)]

	if pairwiseSampleSize > len(aligns):
		logging.warning("The number of available taxonomic \
			matches is smaller than subSampleSize.  \
			Ignoring pairwiseSampleSize.")

	# subsample for pairwise alignments
	subSampledAligns = randomSubSample(aligns, pairwiseSampleSize)

	# convert both pairwise subsample and multiseq iterators to lists
	otuSamples = [str(x.seq) for x in aligns]
	pairwiseSubSamples = [str(x.seq) for x in subSampledAligns]

	pairwiseDist = computePairwiseDistStats(pairwiseSubSamples, distFn)
	otuDist = computeDist(otuSamples, distFn)

	logging.debug("Pairwise Data:\n")
	for key in pairwiseDist.keys():
		logging.debug("%s: %f\n" % (key, pairwiseDist[key]))
	logging.debug("OTU distance: %f\n" % otuDist)
	return(tax, otuDist, pairwiseDist["min"], pairwiseDist["max"], 
		pairwiseDist["avg"])



'''===Notes===
 Directory should look like
 #-/
 #-greg.py
 #-generateTree.py
 #-data/
 #	-Unique_taxonomy.lines
 #	-seq_lin.mapping
 #	-bold_coi_11_05_2015.fasta



"Loading Tree..."
t=loadTree("data/Unique_taxonomy.lines")
"Loading Sequences..."
s=loadSequences("data/seq_lin.mapping")
"Loading BOLD..."
dbFile="data/bold_coi_11_05_2015.fasta"

refDB = indexFasta(dbFile)
distFn="OUTTER_GAP_CONSERVED"

tax="Animalia;Arthropoda;Remipedia;Nectiopoda"
tax="Fungi;Basidiomycota;Agaricomycetes;Boletales"
tax="Protista;Heterokontophyta"
tax="Animalia;Arthropoda;Insecta;Hymenoptera;Eulophidae;Entedoninae;Horismenus;Horismenus opsiphanisDHJ02"
upToLvl=8
"Running query for %s Lvl=%s"  % (tax, upToLvl)
print(computeAlnDistance(t, s, tax, upToLvl, refDB, distFn, "tmp", 15, 10))
'''

