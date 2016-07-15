from Bio import SeqIO, AlignIO
from Bio.Align.Applications import ClustalwCommandline
from generateTree import loadTree, loadSequences, getSequencesOf
from itertools import chain, izip
import re
import subprocess


# filePath: "" path to .aln file
# return: 
def parseALN(filePath):
	handle = open(filePath, "rU")
	records = SeqIO.parse(handle, "clustal")
	handle.close()
	return records



# filePath: "" path to fasta file
# return {id:seq}
def parseFASTA(filePath, nameFn):
	records = SeqIO.index(filePath, "fasta")
	return records



# fastaFilePath: "" path to fasta file
# targetIDs: 	 [] list of IDs to pull
#
# returns: 	 list of SeqIO records
def searchFastaByID(fastaFilePath, targetIDs):
	records = []
	db = parseFASTA(fastaFilePath)
	for record in targetIDs:
	    if record in db.keys():
		records.append(db[record])
	return records



# fastaFile: ""	The fastafile to align.
# takes in fasta file, does alignment, and outputs fastFile.aln
def clustalAlign(fastaFile):
	print("Calling clustalW...")
	outFileName = fastaFile + ".aln"
	subprocess.call(["clustalw", fastaFile, "-outfile=" + outFileName ,"-quiet"])


# line: "" 	A sequence to edit
def replaceInnerGap(line):
	firstLetter = re.search(r'^-*[^-]', line).start()
	lastLetter = re.search(r'[^-]$', line).end()
	line = line[ :firstLetter] + line[firstLetter:lastLetter].replace('-', ':') + line[lastLetter: ]
	return line


# aligns [] 	A list of aligned sequences
# returns 	The mismatch % between aligned sequences
def computeDist_outterGapConserved(aligns):
	miss = 0
	total = len(aligns[0])
	for line in aligns:
		replaceInnerGap(line)
	for col in izip(*aligns):
		colSet = set(col)
		setSize = len(colSet)
		if (setSize > 2) or (setSize == 2 and '-' in colSet):
			miss = miss + 1
	return miss/total


# alginFile: "" An .aln file from clustalw
# distFn: "" 	Uppercase ID for the dist function to use
# returns 	The distance as a % between sequences in an alignFile
def computeDist(alignFile, distFn):
	distFns = {"OUTTER_GAP_CONSERVED":computeDist_outterGapConserved}
	if distFn in distFns.keys():
		return distFns[distFn](alignFile)
	else:
		raise NameError("Error: Invalid grading Fn.  Valid options are:" + distFns.keys())

# t {{}{}} 	A tree from loadTree()
# s {{}{}} 	A tree form loadSequences()
# tax a;b;c;d 	A semicolon-delimited taxonomy string
# upToLvl #	A digit between 1 and 8 indicating depth of tax
# refDB "" 	File path to the fasta DB
# outFile ""	Output name for FASTA of matching Taxonomic entries and their sequences
# distFn ""	Uppercase ID for the dist function to use
def computeConservPct(t, s, tax, upToLvl, refDB, outFile, distFn):
	child_seqs = getSequencesOf(t, s, tax, upToLvl)
	print (child_seqs)
	print("Searching Fasta...")
	children = searchFastaByID(refDB, chain.from_iterable(child_seqs))
	SeqIO.write(children, outFile, "fasta")
	clustalAlign(outFile)
	pct = computeDist(outFile, distFn)
	return pct


#Load Tree and Sequences

t=loadTree("data/Unique_taxonomy.lines")
s=loadSequences("data/seq_lin.mapping")
fastaFilePath="data/bold_coi_11_05_2015.fasta"
outFile="temp.fasta"
tax="Animalia;Arthropoda;Remipedia;Nectiopoda"
upToLvl=8
computeConservPct(t, s, tax, upToLvl, fastaFilePath, outFile, "OUTTER_GAP_CONSERVED")









