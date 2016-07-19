from Bio import SeqIO, AlignIO
from Bio.Align.Applications import ClustalwCommandline, MuscleCommandline
from generateTree import loadTree, loadSequences, getSequencesOf, getLineagesOf
from itertools import chain, izip, combinations

import logging
import numpy as np
import os
import random
import re
import subprocess
import sys

def ensureDirExists(dirPath):
	'''Ensures that a relative directory specified by dirPath exists,\
		by creating it if necessary.

		dirPath "" Relative path to a directory to ensure.
	'''
	# check if temp directorty exists, and create it if not
	tempDir = "temp"
	if not os.path.isdir(tempDir):
		os.mkdir(tempDir)

def parseAln (filePath, alnType="clustal"):
	'''Parses a clustalW .aln file and returns the SeqIO iterator.

	 	filePath  	"" Filepath to .aln file.

	 	return	  	"" A SeqIO iterator.
	'''
	records = SeqIO.parse(open(filePath, "rU"), alnType)
	return records


def indexFasta (filePath, nameFn=lambda x : x.split('|')[0]):
	'''Indexes a fasta file by id, returning a SeqIO index dictionary.

	 	filePath: 	"" Filepath to fasta file.

	 	return 		{SeqIO} A SeqIO index object (non iterable).
	'''
	records = SeqIO.index(filePath, "fasta", key_function=nameFn)
	return records


def searchFastaByID (refDB, targetIDs):
	'''Pulls a list of targetsIDs from a SeqIO index dictionary. 

	 	refDB:  	{seqId: SeqRecord, } SeqIO index object result  
					of indexFasta().
	 	targetIDs: 	["", ] list of IDs to pull.

	 	return	 	[SeqIO.records, ] record.id and record.seq.
	'''
	records = []
	for record in targetIDs:
		# Find the sequence if it exists
		try:
			records.append(refDB[record])
		# Don't care if its not in the DataBase	
		except:
			pass
	return records



def randomSubSample(pool, n):
	'''Randomly subsamples a list without replacement.

	 	n:	# The max number of samples to take.
	 	pool:	[] The pool to take from.
	 
	 	return		A list of up to n items from pool.
	'''
	return random.sample(pool, n)


def muscleAlign(fastaFile):
	'''Takes in fasta file, cals muscle for multi-sequence alignmnet, 
		and outputs "{fastFile}.aln".

	 	fastaFile: 	"" The Filepath to the fastafile to align.
	 
	 	return 		"" The Filepath of the resulting alignment 
					file (fasta format).
	'''

	with open(os.devnull, 'w') as fnull:
		cline = MuscleCommandline(input=fastaFile, out='%s.aln' % fastaFile, clwstrict=True,)
		cline(stderr=fnull, stdout=fnull)
	return "%s.aln" % fastaFile


def clustalAlign(fastaFile):
	'''Takes in fasta file, cals clustalW for multi-sequence alignmnet, 
		and outputs "{fastFile}.aln".

	 	fastaFile: 	"" The Filepath to the fastafile to align.
	 
	 	return 		"" The Filepath of the resulting alignment 
					file (clustW's .aln format).
	'''

	with open(os.devnull, 'w') as fnull:
		clustalWCmd = ['clustalw', '-ALIGN','-INFILE=%s' % fastaFile,
				'-OUTFILE=%s.aln' % fastaFile,
				'-OUTORDER=INPUT']
		subprocess.call(str(clustalWCmd), stderr=fnull, stdout=fnull)
	return "%s.aln" % fastaFile


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
	# The temp directory to write to
	tempDir = "temp"
	ensureDirExists(tempDir)

	# A safe file name delimiter to use as a replacement for unsafe delims.
	fileNameDelim = "_"

	hitsFastaFilePath = "%s/%s.fasta" % (tempDir, 
						tax.replace(';',fileNameDelim))
	child_seqs = chain.from_iterable(getSequencesOf(t, s, tax, upToLvl))
	logging.info("Fetching sequences from database...")
	descendants = list(searchFastaByID(refDB, child_seqs))
	if maxPoolSize > 0:
		if maxPoolSize < len(descendants):
			descendants = randomSubSample(descendants, 
							maxPoolSize)
		else:
			logging.warning(": The number of available taxonomic\
				matches is smaller than maxPoolSize.  \
				Ignoring maxPoolSize.")
	else:
		emptySampleError("maxPoolSize")
		sys.exit()
	# publish matching seqs
	logging.info("Writing results to %s ..." % hitsFastaFilePath)
	SeqIO.write(descendants, open(hitsFastaFilePath, 'w'), "fasta")

	# align seqs
	alignFile = muscleAlign(hitsFastaFilePath)
	aligns = [x for x in parseAln(alignFile)]

	# subsample for pairwise alignments
	if pairwiseSampleSize > 0:
		if pairwiseSampleSize <= len(aligns):
			subSampledAligns = randomSubSample(aligns,
							pairwiseSampleSize)
		else:
			logging.warning("The number of available taxonomic \
				matches is smaller than subSampleSize.  \
				Ignoring pairwiseSampleSize.")
	else:
		emptySampleError("subSampleSize")
		os.exit()

	otuSamples = [str(x.seq) for x in aligns]
	pairwiseSubSamples = [str(x.seq) for x in subSampledAligns]

	logging.info("Computing distance...")
	pairwiseData = computePairwiseDistStats(pairwiseSubSamples, distFn)
	otuDist = computeDist(otuSamples, distFn)

	logging.debug("Pairwise Data:\n")
	for key in pairwiseData.keys():
		logging.debug("%s: %f\n" % (key, pairwiseData[key]))
	logging.debug("OTU distance: %f\n" % otuDist)
	return(tax, otuDist, pairwiseData["min"], pairwiseData["max"], 
		pairwiseData["avg"])

def computeAlnDistAcrossTree(t, maxLvl, minLvl=1):
	total_count=0
	curr_count=0
	for key in t.keys():
		lineages = []
		currentLvl = maxLvl
		while currentLvl >= minLvl: 
			print("============================")
			print("========= Level: %d =========" % (currentLvl) )
			print("============================")
			if currentLvl == 1:
				lineages = [key]
			else:
				lineages = getLineagesOf(t, key, currentLvl)
			i = 1
			for tax in lineages:
				print("%s\t%d/%d\t%d/%d" % \
					(tax, i, len(lineages),
					curr_count, total_count))
				i = i + 1
				curr_count = curr_count+1
			
			currentLvl = currentLvl - 1

	
'''===Notes===
 Directory should look like
 #-/
 #-greg.py
 #-generateTree.py
 #-data/
 #	-Unique_taxonomy.lines
 #	-seq_lin.mapping
 #	-bold_coi_11_05_2015.fasta

'''
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
upToLvl=8

"Running query for %s Lvl=%s"  % (tax, upToLvl)
print(computeAlnDistance(t, s, tax, upToLvl, refDB, distFn, "tmp", 15, 10))

