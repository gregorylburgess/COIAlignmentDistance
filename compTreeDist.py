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
            samples = []
            samplePool = []
            rslts = []
            distFns = []
            outDirs = []
            if currentLvl == 1:
                lineages = [key]
            else:
                lineages = getLineagesOf(t, key, currentLvl)
            "%s %d" % (key, currentLvl)
            distFns = [distFn for lineage in lineages]
            outDirs = [outDir for lineage in lineages]
            for tax in lineages:
                try:
                    samplePool = subsampleTaxa(t, s, tax, 8, refDB, distFn, outDir)

                except:
                    samplePool = []
                    err.write("%s\t%s\n" % (tax, sys.exc_info()))
                samples.append(samplePool)
            try:
                print("Starting Threads")
                print(lineages)
                pool = Pool(processes=poolSize)
                rslts = pool.map(unwrapAndCall,
                        [(computeAlnDistance, tax, samplePool, distFn, outDir) \
                         for tax, samplePool, distFn, outDir in zip(lineages, samples, distFns, outDirs)])
                print("Collecting Threads")
                pool.close()
                pool.join()
                print("Done Collecting Threads")
            except:
                print(sys.exc_info())


            i=0
            rsltLength = len(rslts)
            for rslt in rslts:
                i = i + 1
                print("%d/%d" % (i, rsltLength))
                try:
                    out.write("%s\t%f\t%f\t%f\t%f\n" % rslt)
                except:
                    print(sys.exc_info())
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
    print "Starting with %d threads..."%poolSize
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
                    print "%s %d" % (key, currentLvl)
                    lineages = getLineagesOf(t, key, currentLvl)

                outfiles = ["%s/%s.fasta" % (outDir, sanitizeFileName(tax, "_")) for tax in lineages]
                child_seqs = [chain.from_iterable(getSequencesOf(t, s, tax)) for tax in lineages]
                for seq in child_seqs:
                    print (seq)
                descendants = [list(searchFastaByID(refDB, seqs)) for seqs in child_seqs]
                if poolSize > 1:
                    try:
                        rslt = pool.map(unwrapAndCall, [(alignDescendants, seqs, outfile) for seqs, outfile in zip(descendants, outfiles)])

                    except:
                        print(sys.exc_info())
                else:
                    for tax in lineages:
                        try:
                            print tax
                            alignDescendants(t, s, refDB, tax, outDir, maxLvl)
                            out.write("%s\n" % tax)
                        except:
                            err.write("%s\t%s\n" % (tax, sys.exc_info()))

    except KeyboardInterrupt:
        pool.terminate()
        return

    "Done!"
    err.close()
    out.close()



# computeAlnDistAcrossTree(8, 1, poolSize=4)
if len(sys.argv) < 4:
    "Usage: (align|dist) #threads outDir"
if sys.argv[1] == "align":
    alignDescendantsAcrossTree(7, 0, int(sys.argv[2]), sys.argv[3])
elif sys.argv[1] == "dist":
    computeAlnDistAcrossTree(7, 0, int(sys.argv[2]), sys.argv[3])