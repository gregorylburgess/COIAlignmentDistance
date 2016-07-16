from Bio import SeqIO, AlignIO
from Bio.Align.Applications import ClustalwCommandline
from generateTree import loadTree, loadSequences, getSequencesOf
from itertools import chain, izip, combinations
import os
import random
import re
import subprocess

# TODO
# Check that the subsample is <= the number of sequences in that node



# Parses a clustalW .aln file and returns the SeqIO iterator.
# 	filePath  	"" Filepath to .aln file.
#
# 	return	  	"" A SeqIO iterator.
def parseAln (filePath):
	handle = open(filePath, "rU")
	records = SeqIO.parse(handle, "clustal")
	return records


# Indexes a fasta file by id, returning a SeqIO index dictionary.
# 	filePath: 	"" Filepath to fasta file.
#
# 	return 		{SeqIO} A SeqIO index object (non iterable).
def indexFasta (filePath, nameFn=lambda x : x.split('|')[0]):
	records = SeqIO.index(filePath, "fasta", key_function=nameFn)
	return records


# Pulls a list of targetsIDs from a SeqIO index dictionary. 
# 	refDB:  	 {SeqIO} SeqIO index object result of indexFasta().
# 	targetIDs: 	 [""] list of IDs to pull.
#
# 	return	 	 [SeqIO.records] record.id and record.seq.
def searchFastaByID (refDB, targetIDs):
	records = []
	for record in targetIDs:
		# Find the sequence if it exists
		try:
			records.append(refDB[record])
		# Don't care if its not in the DataBase	
		except:
			pass
	return records


# Randomly subsamples a list without replacement
#	n:	# The max number of samples to take.
#	pool:	[] The pool to take from
#
#	return		A list of up to n items from pool.
def randomSubSample(pool, n):
	return random.sample(pool, n)


# Takes in fasta file, cals clustalW for multi-sequence alignmnet, and outputs "{fastFile}.aln"
# 	fastaFile: 	"" The Filepath to the fastafile to align.
#
# 	return 		"" The Filepath of the resulting alignment file (clustW's .aln format)
def clustalAlign (fastaFile):
	print("Running clustalW...")
	outFileName = fastaFile + ".aln"
	subprocess.call(["clustalw", fastaFile, "-outfile=" + outFileName ,"-quiet"])
	return outFileName

# Replaces the internal gaps '-' of a string by ':'
# 	line:  	"" A sequence to edit
#
# 	return 	"" The string with inner gaps replaced by ':'
def replaceInnerGap (line):
	frstA = re.search(r'^-*[^-]', line).end() - 1
	lastA = re.search(r'[^-]-*$', line).start()
	return line[ :frstA] + (line[frstA:lastA].replace('-', ':')) + line[lastA: ]

# Computes the distance between the sequences in a .aln file, counting padding gaps as wildcard matches.
# 	aligns	 	[""] List of seqs to align"" File path to the .aln file from clustalw
#
# 	returns 	The mismatch % between aligned sequences
def computeDist_outterGapConserved (aligns):
	miss = 0
	succ = 0
	for i in range(len(aligns)):
		aligns[i] = replaceInnerGap(aligns[i]) 

	for col in izip(*aligns):
		colSet = set(col)
		setSize = len(colSet)
		if ':' in colSet or (setSize > 2) or (setSize == 2 and not '-' in colSet):
			miss = miss + 1
		else:
			succ = succ + 1
	#print ("Miss=" + str(miss))
	#print ("Hit="+str(succ))
	return miss / (miss + float(succ))

# Computes the distance between sequences in a .aln file, using the distFn specified.
# 	algins: 	[] List of seqs to align"" File path to the .aln file from clustalw
# 	distFn: 	"" Uppercase ID for the dist function to use
#	
# 	return		The distance as a % between sequences in an alignFile
def computeDist (aligns, distFn):
	distFns = {"OUTTER_GAP_CONSERVED": computeDist_outterGapConserved}
	
	if distFn in distFns.keys():
		return distFns[distFn](aligns)
	else:
		raise NameError("Error: Invalid grading Fn.  Valid options are:" + distFns.keys())

# Does pairwise subsampling on a list of sequences.  Returns a list of sample diversity (as % difference) rates.
#	aligns:		[] The sequences to pairwise align.
def pairwiseSubSample(alignedSeqs, distFn):
	sampleSize = len(alignedSeqs)
	dists = [] 
	for i,j in combinations(range(sampleSize),2):
			dists.append( computeDist([str(alignedSeqs[i]), str(alignedSeqs[j])], distFn))
	return dists

# Computes the min, max, avg % genetic variance amongst members of a given OTU.  Primary function.
# 	t		{{}{}} A tree from loadTree()
# 	s		{{}{}} A tree form loadSequences()
# 	refDB 		{} A SeqIO fast index.  Result of calling indexFasta()
# 	tax		"a;b;c;d" A semicolon-delimited taxonomy string
# 	upToLvl		# A digit between 1 and 8 indicating depth of tax
# 	distFn		"" Uppercase ID for the dist function to use
#	maxPoolSize	# A number indicating the maxmimum size of a the taxa to multi-seq align.
#			Default = no sub-sampling.
#	pairwiseSampleSize	# A number indicating the maxmimum size of a subsample from the tax pool to pariwise-seq align.
#			Default = no sub-sampling.
# 	returns		#.## Percentage diversity between sequences at a given taxanomic level
def computeOtuDistance (t, s, tax, upToLvl, refDB, distFn, maxPoolSize = 0, subSampleSize=0):
	hitsFastaFilePath = "data/" + tax + ".fasta"
	child_seqs = chain.from_iterable (getSequencesOf (t, s, tax, upToLvl))
	print("Fetching sequences from database...")
	children = list(searchFastaByID (refDB, child_seqs))
	if maxPoolSize > 0:
		children = randomSubSample (children, maxPoolSize)

	# publish matching seqs
	print("Writing results to " + hitsFastaFilePath)
	SeqIO.write(children, hitsFastaFilePath, "fasta")

	# align seqs
	alignFile = clustalAlign(hitsFastaFilePath)
	aligns = [x for x in parseAln(alignFile)]

	# subsample for pairwise alignments
	if subSampleSize > 0:
		alignsSubSample = randomSubSample (aligns, subSampleSize)

	cleanedAligns = [str(x.seq) for x in alignsSubSample]
	#=====================works ok up to here==============
	# at this point cleanedAligns is a cleaned list of strings representing the sequences that are ready to be pairwise aligned

	#compute distances
	print("Computing distance.")
	dists = pairwiseSubSample(alignsSubSample, distFn)
	print(dists)
	max_value = max(dists)
	min_value = min(dists)
	avg_value = sum(dists)/len(dists)
	print(min_value, avg_value, max_value)

#directory should look like
#-/
#-greg.py
#-generateTree.py
#-data/
#	-Unique_taxonomy.lines
#	-seq_lin.mapping
#	-bold_coi_11_05_2015.fasta

#Load Tree and Sequences


# print("Loading Tree...")
# t=loadTree("data/Unique_taxonomy.lines")
# print("Loading Sequences...")
# s=loadSequences("data/seq_lin.mapping")
# print("Loading BOLD...")
# dbFile="data/bold_coi_11_05_2015.fasta"

# refDB = indexFasta(dbFile)
# distFn="OUTTER_GAP_CONSERVED"

# tax="Animalia;Arthropoda;Remipedia;Nectiopoda"
# upToLvl=8

# print "Running query for %s Lvl=%s"  % ( tax, upToLvl)

# print  computeOtuDistance(t, s, tax, upToLvl, refDB, distFn, 10, 3)

