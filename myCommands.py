import os
from Bio import SeqIO, SearchIO
import sys
import glob
from Bio import SeqRecord
from Bio.SeqUtils import six_frame_translations
import subprocess
import itertools
import numpy as np
from Bio import AlignIO


# generateProtNucFiles(nucFile, outDir)  generate protein translation file from the fasta file downloaded from BOLD
# runAlignments(myFileProt, outDir) aling the protein seqeunces and generate the pal2nal file
# computeDist(seq1, seq2) compute effective distance between two sequences seq1 and seq2
# computeNucAliStats(outPal2NalFile) computes nucl. stats using pal2nal file (myMin, myMax, myAvg, myStd)
# validCols, probs = runHMMer(proteinAliFName) returns the valid HMM columns and the probabilities dict for the valid columns
# computeHmmStats(aliFName, validCols) returns a dict of all HMMer bitscore for the sequences


def generateProtNucFiles(fileName, outDir):
    reads = SeqIO.parse(fileName, 'fasta')
    prots = []
    dna = []
    for read in reads:
        # there are sequence from markers other than COI. Filter those out
        # if read.description.find("COI-5P") == -1:
        #     print "ignoring %s " % read.description
        #     continue
        sequence = read.seq.ungap("-")
        translation = sequence.translate(table=9)
        if translation.count("*"):
            hasGoodTr = False
            # find other possible reading frames or die trying
            for frame in [1,2]:
                s = sequence[frame:]
                t = s.translate(table=9)
                if not t.count("*"):
                    sequence = s
                    translation = t
                    hasGoodTr = True
                    break
            if not hasGoodTr:
                print "\n\n\n\n"
                print "****** seqeuence %s contains a stop in its translation" % read.id
                print "\n\n\n\n"
                print translation
                sys.exit(0)
        dna.append(SeqRecord.SeqRecord(id=read.id, seq=sequence, name=""))
        prots.append(SeqRecord.SeqRecord(id=read.id, seq=translation, name=""))
    print "Writing %s sequences for file %s" % (len(prots), fileName)
    SeqIO.write(dna, open(os.path.join(outDir,os.path.splitext(os.path.basename(fileName))[0]+".fna"), 'w'), 'fasta')
    SeqIO.write(prots, open(os.path.join(outDir,os.path.splitext(os.path.basename(fileName))[0]+".faa"), 'w'), 'fasta')

def computeDist(seq1, seq2):
    import sys
    if len(seq1) != len(seq2):
        # TODO CHANGE TO EXCEPTION
        print "Seqeunces are not of the same length"
        sys.exit(0)
    start = 0
    end = len(seq1)

    # find true start
    for x in range(len(seq1)):
        if seq1[x] != '-':
            start = x
            break
    for x in range(len(seq2)):
        if seq2[x] != '-':
            if start < x:
                start = x 
            break
    for x in range(len(seq1))[::-1]:
        if seq1[x] != '-':
            end = x
            break
    for x in range(len(seq2))[::-1]:
        if seq2[x] != '-':
            if start > x:
                end = x 
            break
    dist = sum([1 for x in range(start, end+1) if seq1[x] != seq2[x] ])
    print seq1.id, seq2.id, start, end, dist
    return dist
        
def computeNucAliStats(outPal2NalFile):
    # TODO: move this top
    ali = AlignIO.read(outPal2NalFile, 'clustal')
    distMatrix = np.zeros((len(ali), len(ali)))

    for pair in itertools.combinations(range(len(ali)), 2):
        distMatrix[pair[0], pair[1]] = computeDist(ali[pair[0]], ali[pair[1]])

    print distMatrix
    myAvg = np.average(distMatrix)
    myStd = np.std(distMatrix)
    myMax = distMatrix[np.unravel_index(distMatrix.argmax(), distMatrix.shape)]
    # fill diagonal with max to ignore z
    np.fill_diagonal(distMatrix, distMatrix.max())
    myMin = distMatrix[np.unravel_index(distMatrix.argmin(), distMatrix.shape)]
    return (myMin, myMax, myAvg, myStd)
        
def runAlignments(myFileProt, outDir):
    # assumes that nucleotide and amino acid files have the smae name scheme except for extension
    # faa versus fna
    
    baseName = os.path.splitext(os.path.basename(myFileProt))[0]

    myFileDna = os.path.join(os.path.dirname(myFileProt), baseName+".fna")
    protAliOutFile = os.path.join(outDir, baseName+".aln")
    pal2nalOutFile = os.path.join(outDir, baseName+".pal2nal")
    
    sequences = SeqIO.to_dict(SeqIO.parse(myFileProt, 'fasta'))
    if len(sequences) == 1:
        # nothing to align, just write to files and continute to next file
        SeqIO.write(sequences.values(), open( protAliOutFile, 'w'), 'clustal')
        sequences = SeqIO.to_dict(SeqIO.parse(myFileDna, 'fasta'))
        SeqIO.write(sequences.values(), open(pal2nalOutFile, 'w'), 'clustal')
        print "Single seq file, exiting"
        return
    clustalWCmd = ['/Users/mahdi/programs/clustalw-2.1-macosx/clustalw2', '-ALIGN',
                   '-INFILE=%s' % myFileProt, '-TYPE=PROTEIN', '-OUTFILE=%s' % protAliOutFile, '-OUTORDER=INPUT']
    # mafftCmd = ['mafft-linsi', '--amino', '--clustalout', myFileProt , '>', protAliOutFile]


    pal2nalcmd = ['perl', '/Users/mahdi/programs/pal2nal.v14/pal2nal.pl', protAliOutFile, myFileDna, '-codontable', '9',  ]
    with open(os.devnull, 'w') as fnull, open(pal2nalOutFile, 'w') as out:
        subprocess.call(clustalWCmd, stderr=fnull, stdout=fnull)
        subprocess.call(pal2nalcmd, stdout=out)
    print "completed file %s" % baseName
    
def runHMMer(proteinAliFName):

    hmmerCmd = ["/Users/mahdi/programs/hmmer-3.1b2-macosx-intel/binaries/hmmbuild",  
    "-O", '/tmp/effectiveAli',  "%s.hmm" % proteinAliFName,  "%s" % proteinAliFName]
    hmmLogo = ["/Users/mahdi/programs/hmmer-3.1b2-macosx-intel/binaries/hmmlogo", "%s.hmm" % proteinAliFName]

    aminoAcids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
    print hmmerCmd
    with open(os.devnull, 'w') as fnull, open("%s.logo" %  proteinAliFName, 'w') as logoFile:
        subprocess.call(hmmerCmd, stderr=fnull, stdout=fnull)
        subprocess.call(hmmLogo, stderr=fnull, stdout=logoFile)

    concensusLine  = ""
    for line in open('/tmp/effectiveAli','r'):
        line = line.rstrip()
        if line[0:4] == "#=GC":
            concensusLine += line.split()[2]
    effectiveCols = [x for x in range(len(concensusLine)) if concensusLine[x]=='x']
    logoFile = open("%s.logo" % proteinAliFName, 'r')
    logoFile.readline()
    logoFile.readline()

    probs ={}
    for line in logoFile:
        line = line.rstrip()

        data = line.split()

        colNum = data[0].replace(":", "")
        if data[0][0] not in ['0','1','2','3','4','5','6','7','8','9']:
            break
        total = sum([float(x) for x in data[1:-2]])
        print line
        normalizedPro = [float(x)/total for x in data[1:-2]]
        probs[int(colNum)] = dict(zip(aminoAcids, normalizedPro))
    return effectiveCols, probs

def computeStats(aliFName, validCols):

    # Quick and VERY slow
    # convert o matrix and operate directly on it
    ali = AlignIO.read(aliFName, 'clustal')
    scores = {}

    for seqNumRemove in range(len(ali)):
        newSeqs = AlignIO.MultipleSeqAlignment([])
        aliWithValidCols = AlignIO.MultipleSeqAlignment([])

        for seqNum in range(len(ali)):
            if seqNum != seqNumRemove:
                newSeqs.append(ali[seqNum])

        aliWithValidCols = newSeqs[:,validCols[0]:validCols[0]+1]
        querySeq = ali[seqNumRemove, validCols[0]:validCols[0]+1]

        for col in validCols[1:]:
            aliWithValidCols += newSeqs[:,col:col+1]
            querySeq += ali[seqNumRemove, col:col+1]
        AlignIO.write(aliWithValidCols, open("/tmp/tempPartialAli", 'w'), 'clustal')
        SeqIO.write(querySeq, open('/tmp/querySeq', 'w'), 'fasta')
        print "seq %s out of %s completed" % (seqNumRemove, len(ali))
        print "running hmmer and gathering stats"

        hmmBuildCmd = ["/Users/mahdi/programs/hmmer-3.1b2-macosx-intel/binaries/hmmbuild", "/tmp/tempPartialAli.hmm", "/tmp/tempPartialAli"]
        hmmScanCmd  = ["/Users/mahdi/programs/hmmer-3.1b2-macosx-intel/binaries/hmmscan",
                      "/tmp/tempPartialAli.hmm", "/tmp/querySeq"]
        hmmPressCmd = ["/Users/mahdi/programs/hmmer-3.1b2-macosx-intel/binaries/hmmpress",
                      "/tmp/tempPartialAli.hmm"]
        with open(os.devnull, 'w') as fnull, open("/tmp/alResult", 'w') as out:
            subprocess.call(hmmBuildCmd, stdout=fnull)
            subprocess.call(hmmPressCmd, stderr=fnull, stdout=fnull)

            subprocess.call(hmmScanCmd, stderr=fnull, stdout=out)


            print hmmBuildCmd
            print hmmScanCmd
        search =  SearchIO.parse("/tmp/alResult", "hmmer3-text").next()
        bScore = search.hsps[0].bitscore
        scores[search.id] = bScore
    return scores
