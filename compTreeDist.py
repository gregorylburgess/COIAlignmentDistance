from Bio import SeqIO, AlignIO
from Bio.Align.Applications import ClustalwCommandline, MuscleCommandline
from generateTree import loadTree, loadSequences, getSequencesOf, getLineagesOf
from helpers import *
from computeDist import *
from multiprocessing import Pool
import sys


def cleanUp(outDir, taxFileName):
	for fileExt in ["fasta","fasta.aln"]:
		fileName = "%s/%s.%s" % (outDir, taxFileName, fileExt)
		if os.path.isfile(fileName):
			os.remove(fileName)

def computeAlnDistAcrossTree(maxLvl=8, minLvl=0, poolSize=4):
	treeFile = "data/Unique_taxonomy.lines"
	seqsFile = "data/seq_lin.mapping"
	dbFile = "data/bold_coi_11_05_2015.fasta"
	distFn = "OUTTER_GAP_CONSERVED"
	outDir = "tree"
	maxPoolSize = 10
	pairwiseSampleSize = 5
	"Loading..."
	t = loadTree(treeFile)
	s = loadSequences(seqsFile)
	refDB = indexFasta(dbFile)
	ensureDirExists(outDir)
	pool = Pool(processes=poolSize)
	f = open("%s/%s.txt" % (outDir, "metaData"),'w')
	f.write("treeFile:\t%s\n\
		seqsFile:\t%s\n\
		dbFile:\t%s\n\
		distFn:\t%s\n\
		outDir:\t%s\n\
		maxPoolSize:\t%s\n\
		pairwiseSampleSize:\t%s\n" %  (treeFile, seqsFile, dbFile, 
			distFn, outDir, maxPoolSize, pairwiseSampleSize))
	f.close() 
	
	err = open("%s/ERRORS.txt" % (outDir),'w')
	err.write("")
	err.close() 
	err = open("%s/ERRORS.txt" % (outDir),'a')
	"Stating..."
	out = open("%s/rslt.txt" % outDir,'w')
	for key in t.keys():
		lineages = []
		currentLvl = minLvl
		while currentLvl <= maxLvl: 
			currentLvl = currentLvl + 1
			if currentLvl == 1:
				lineages = [key]
			else:
				"%s %d" % (key, currentLvl)
				lineages = getLineagesOf(t, key, currentLvl)
			for tax in lineages:
				try:
					rslt = computeAlnDistance(t, s, tax, 8, refDB, 
						distFn, outDir, maxPoolSize, 
						pairwiseSampleSize)
					cleanUp(outDir, sanitizeFileName(tax,"_"))
					out.write("%s\t%f\t%f\t%f\t%f\n" % rslt)
				except:
						err.write("%s\t%s\n" % (tax,sys.exc_info()[0]))
						
			#args = [(t, s, tax, 8, refDB, distFn, 
			#	outDir, maxPoolSize, pairwiseSampleSize) for tax in lineages]
			#pool.map_async(compute, args)
	"Done!"
	err.close() 
	out.close()

computeAlnDistAcrossTree(8, 1, poolSize=4)
