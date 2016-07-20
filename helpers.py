import os
import random
from Bio import SeqIO, AlignIO
from Bio.SeqIO import *
from Bio.Align.Applications import ClustalwCommandline, MuscleCommandline
from generateTree import loadTree, loadSequences, getSequencesOf, getLineagesOf


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


def sanitizeFileName(fileName, replacement):
	fileName = fileName.replace(' ',replacement)
	fileName = fileName.replace(';',replacement)
	fileName = fileName.replace('/',replacement)
	fileName = fileName.replace('.',replacement)
	return fileName

def randomSubSample(pool, n):
	'''Randomly subsamples a list without replacement.

	 	n:	# The max number of samples to take.
	 	pool:	[] The pool to take from.
	 
	 	return		A list of up to n items from pool.
	'''
	k = min(len(pool), n)
	return random.sample(pool, k)


def muscleAlign(fastaFile):
	'''Takes in fasta file, cals muscle for multi-sequence alignmnet, 
		and outputs "{fastFile}.aln".

	 	fastaFile: 	"" The Filepath to the fastafile to align.
	 
	 	return 		"" The Filepath of the resulting alignment 
					file (fasta format).
	'''

	with open(os.devnull, 'w') as fnull:
		cline = MuscleCommandline(input=fastaFile, out='%s.aln' % \
					fastaFile, clwstrict=True)
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

