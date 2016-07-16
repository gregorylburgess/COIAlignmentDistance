from Bio import SeqIO, AlignIO
from Bio.Align.Applications import ClustalwCommandline
from generateTree import loadTree, loadSequences, getSequencesOf
from itertools import chain, izip, combinations

# TOD: this only exists in version 3.4
# preferably work with numpy to remove any dependency on python 3
# from statistics import mean, stdev

# TODO: Add a DEBUG boolean and only use print when DEBUG is set to true.


import numpy as np
import os
import random
import re
import subprocess

# TODO: 
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

	 	refDB:  	{seqId: SeqRecord, } SeqIO index object result of 
					indexFasta().
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
		subprocess.call(clustalWCmd, stderr=fnull, stdout=fnull)
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

	 	aligns	 	[""] List of seqs to align"" File path to the 
					.aln file from clustalw.
	 
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

	 	aligns: 	[] List of seqs to align"" File path to the 
					.aln file from clustalw
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

	print distMatrix
	myMin = min(valList)
	myMax = max(valList)
	myAvg = np.mean(valList)
	myStd = np.std(valList)
	sol = (myMin, myMax, myAvg, myStd)
	# TODO: use DEBUG to print
	print 'min: %f\nmax: %f\navg: %f\nstd: %f' % sol
	return sol


def computeAlnDistance (t, s, tax, upToLvl, refDB, distFn, maxPoolSize, \
				pairwiseSampleSize):
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
	 	maxPoolSize	# A number indicating the maxmimum size of a
					 the taxa to multi-seq align.
	 			Default = no sub-sampling.
	 	pairwiseSampleSize	# A number indicating the maxmimum 
						size of a subsample from the 
						tax pool to pariwise-seq align.
	 			Default = no sub-sampling.
	 	returns		#.## Percentage diversity between sequences 
					at a given taxanomic level
	'''
	# TODO: I don't think this conditions is needed.
	if pairwiseSampleSize > maxPoolSize:
		raise Exception("subSampleSize must be smaller than\
					 maxPoolSize.")

	
	# TODO: create a temp directory specifically for results to avoid polluting the
	# data directory
	hitsFastaFilePath = "data/%s.fasta" % tax
	child_seqs = chain.from_iterable(getSequencesOf(t, s, tax, upToLvl))
	# TODO: use DEBUG to print
	print("Fetching sequences from database...")
	descendants = list(searchFastaByID(refDB, child_seqs))
	if maxPoolSize > 0:
		if maxPoolSize < len(descendants):
			descendants = randomSubSample(descendants, 
							maxPoolSize)
		else:
			print "WARNING: The number of available taxonomic\
				matches is smaller than maxPoolSize.  \
				Ignoring maxPoolSize."
	else:
		# TODO: should stop otherwise
		pass
	# publish matching seqs
	# TODO: use DEBUG to print
	print "Writing results to %s " % hitsFastaFilePath
	SeqIO.write(descendants, open(hitsFastaFilePath, 'w'), "fasta")

	# align seqs
	alignFile = clustalAlign(hitsFastaFilePath)
	aligns = [x for x in parseAln(alignFile)]

	# subsample for pairwise alignments
	if pairwiseSampleSize > 0:
		if pairwiseSampleSize <= len(aligns):
			subSampledAligns = randomSubSample(aligns,
							pairwiseSampleSize)
		else:
			print "WARNING: The number of available taxonomic \
				matches is smaller than subSampleSize.  \
				Ignoring pairwiseSampleSize."
	else:
		# TODO: Same as above. Nothign to compute if this pairwiseSampleSize <= 0

	cleanedAligns = [str(x.seq) for x in subSampledAligns]
	# TODO: use DEBUG to print
	print "Computing distance."
	return computePairwiseDistStats(cleanedAligns, distFn)

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

# print "Loading Tree..."
# t=loadTree("data/Unique_taxonomy.lines")
# print "Loading Sequences..."
# s=loadSequences("data/seq_lin.mapping")
# print "Loading BOLD..."
# dbFile="data/bold_coi_11_05_2015.fasta"

# refDB = indexFasta(dbFile)
# distFn="OUTTER_GAP_CONSERVED"

# tax="Animalia;Arthropoda;Remipedia;Nectiopoda"
# upToLvl=8

# print "Running query for %s Lvl=%s"  % ( tax, upToLvl)

# print  computeAlnDistance(t, s, tax, upToLvl, refDB, distFn, 15, 10)
