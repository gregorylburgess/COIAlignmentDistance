from Bio import SeqIO, AlignIO
from Bio.Align.Applications import ClustalwCommandline, MuscleCommandline
from generateTree import loadTree, loadSequences, getSequencesOf, getLineagesOf
from helpers import *
from computeDist import *
from multiprocessing import Pool
import sys


def cleanUp(outDir, taxFileName):
    for fileExt in ["fasta", "fasta.aln"]:
        fileName = "%s/%s.%s" % (outDir, taxFileName, fileExt)
        if os.path.isfile(fileName):
            os.remove(fileName)


def computeAlnDistAcrossTree(maxLvl=8, minLvl=0, poolSize=4, outDir="tree"):
    treeFile = "data/Unique_taxonomy.lines"
    seqsFile = "data/seq_lin.mapping"
    dbFile = "data/bold_coi_11_05_2015.fasta"
    distFn = "OUTTER_GAP_CONSERVED"
    maxPoolSize = 10
    pairwiseSampleSize = 5
    "Loading..."
    t = loadTree(treeFile)
    s = loadSequences(seqsFile)
    refDB = indexFasta(dbFile)
    ensureDirExists(outDir)
    pool = Pool(processes=poolSize)
    f = open("%s/%s.txt" % (outDir, "metaData"), 'w')
    f.write("treeFile:\t%s\n\
		seqsFile:\t%s\n\
		dbFile:\t%s\n\
		distFn:\t%s\n\
		outDir:\t%s\n\
		maxPoolSize:\t%s\n\
		pairwiseSampleSize:\t%s\n" % (treeFile, seqsFile, dbFile,
                                      distFn, outDir, maxPoolSize, pairwiseSampleSize))
    f.close()

    err = open("%s/ERRORS.txt" % (outDir), 'w')
    err.write("")
    err.close()
    err = open("%s/ERRORS.txt" % (outDir), 'a')
    "Stating..."
    out = open("%s/rslt.txt" % outDir, 'w')
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
                    cleanUp(outDir, sanitizeFileName(tax, "_"))
                    out.write("%s\t%f\t%f\t%f\t%f\n" % rslt)
                except:
                    err.write("%s\t%s\n" % (tax, sys.exc_info()[0]))

                # args = [(t, s, tax, 8, refDB, distFn,
                #	outDir, maxPoolSize, pairwiseSampleSize) for tax in lineages]
                # pool.map_async(compute, args)
    "Done!"
    err.close()
    out.close()


def write(x, y):
    print(x)
    print(y)


def alignDescendantsAcrossTree(maxLvl=8, minLvl=0, poolSize=4, outDir="tree"):
    treeFile = "data/Unique_taxonomy.lines"
    seqsFile = "data/seq_lin.mapping"
    dbFile = "data/bold_coi_11_05_2015.fasta"
    "Loading..."
    t = loadTree(treeFile)
    s = loadSequences(seqsFile)
    refDB = indexFasta(dbFile)
    ensureDirExists(outDir)
    pool = Pool(processes=poolSize)

    err = open("%s/ERRORS.txt" % (outDir), 'w')
    err.write("")
    err.close()
    err = open("%s/ERRORS.txt" % (outDir), 'a')
    "Starting..."
    out = open("%s/rslt.txt" % outDir, 'w')
    try:
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
                    print(tax)
                    out.write("%s\n" % tax)
                if poolSize > 1:
                    try:
                        pool.imap_unordered(unwrapAndCall,
                                            [(alignDescendants, t, s, refDB, tax, outDir, maxLvl) for x in lineages])

                    except:
                        err.write("%s\t%s\n" % (tax, sys.exc_info()[0]))
                else:
                    for tax in lineages:
                        try:
                            print tax
                            alignDescendants(t, s, refDB, tax, outDir, maxLvl)
                            out.write("%s\n" % tax)
                        except:
                            err.write("%s\t%s\n" % (tax, sys.exc_info()[0]))

    except KeyboardInterrupt:
        pool.close()
        pool.join()
    "Done!"
    err.close()
    out.close()


# computeAlnDistAcrossTree(8, 1, poolSize=4)
if len(sys.argv) < 4:
    "Usage: (align|dist) #threads outDir"
if sys.argv[1] == "align":
    alignDescendantsAcrossTree(8, 0, sys.argv[2], sys.argv[3])
elif sys.argv[1] == "dist":
    computeAlnDistAcrossTree(8, 0, sys.argv[2], sys.argv[3])