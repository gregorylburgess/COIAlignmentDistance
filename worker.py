from Bio import SeqIO, AlignIO
from Bio.Align.Applications import ClustalwCommandline, MuscleCommandline
from generateTree import loadTree, loadSequences, getSequencesOf, getLineagesOf
from helpers import *
from computeDist import *
from multiprocessing import Pool

def compute(t, s, tax, upToLvl, refDB, distFn, outDir, maxPoolSize, 
		pairwiseSampleSize):
	
	print(tax+" "+str(upToLvl))
	sol = computeAlnDistance(t, s, tax, upToLvl, refDB, distFn, outDir, 
				maxPoolSize, pairwiseSampleSize)
	print (sol)
	taxFileName = tax.replace(";","_")
	f = open("%s/%s.txt" % (outDir, taxFileName),'w')
	f.write("%s\t%f\t%f\t%f\t%f" % sol) 
	f.close() 
	cleanUp(outDir, taxFileName)


def cleanUp(outDir, taxFileName):
	for fileExt in ["fasta","fasta.aln"]:
		fileName = "%s/%s.%s" % (outDir, taxFileName, fileExt)
		if os.path.isfile(fileName):
			os.remove(fileName)

def computeAlnDistAcrossTree(maxLvl, minLvl=0, poolSize=4):
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
	"Stating..."
	for key in t.keys():
		lineages = []
		currentLvl = minLvl
		while currentLvl <= maxLvl: 
			currentLvl = currentLvl + 1
			if currentLvl == 1:
				lineages = [key]
			else:
				print(key+" "+str(currentLvl))
				lineages = getLineagesOf(t, key, currentLvl)
			for tax in lineages:
				compute(t, s, tax, currentLvl, refDB, distFn, 
					outDir, maxPoolSize, pairwiseSampleSize)
			#args = [(t, s, tax, 8, refDB, distFn, 
			#	outDir, maxPoolSize, pairwiseSampleSize) for tax in lineages]
			#pool.map_async(compute, args)
	"Done!"
	f = open("%s/%s.txt" % (outDir, "DONE"),'w')
	f.write("DONE")
	f.close() 

