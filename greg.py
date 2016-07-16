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
	return records



# filePath: "" path to fasta file
# return { x.id : x.seq } A SeqIO index object (non iterable)
def indexFASTA(filePath, nameFn=lambda x : x.split('|')[0]):
	records = SeqIO.index(filePath, "fasta", key_function=nameFn)
	return records



# fastaFilePath: "" path to fasta file
# targetIDs: 	 [] list of IDs to pull
#
# returns: 	 [SeqIO.record, ...] record.id and record.seq
def searchFastaByID(fastaFilePath, targetIDs):
	records = []
	db = indexFASTA(fastaFilePath)
	for record in targetIDs:#TODO split record names on |, take [0]
		# Find the sequence if it exists
		try:
			records.append(db[record])
		# Don't care if its not in the DataBase	
		except:
			pass
	return records

 

# fastaFile: ""	The fastafile to align.
# takes in fasta file, does alignment, and outputs fastFile.aln
# returns the filepath of the resulting alignment file (clustW's .aln format)
def clustalAlign(fastaFile):
	print("Calling clustalW...")
	outFileName = fastaFile + ".aln"
	subprocess.call(["clustalw", fastaFile, "-outfile=" + outFileName ,"-quiet"])
	return outFileName


# line: "" 	A sequence to edit
# return ""	The string with inner gaps replaced by ':'
def replaceInnerGap(line):
	frstA = re.search(r'^-*[^-]', line).end() - 1
	lastA = re.search(r'[^-]-*$', line).start()
	return line[ :frstA] + (line[frstA:lastA].replace('-', ':')) + line[lastA: ]


# aligns [] 	A list of aligned sequences
# returns 	The mismatch % between aligned sequences
def computeDist_outterGapConserved(alignFile):
	miss = 0
	succ = 0
	aligns = [replaceInnerGap(str(x.seq)) for x in parseALN(alignFile)]
	for col in izip(*aligns):
		colSet = set(col)
		print(colSet)
		setSize = len(colSet)
		if ':' in colSet or (setSize > 2) or (setSize == 2 and not '-' in colSet):
			print ("M")
			miss = miss + 1
		else:
			succ = succ + 1
	print ("Miss=" + str(miss))
	print("Hit="+str(succ))
	return miss / (miss + float(succ))


# alginFile: "" An .aln file from clustalw
# distFn: "" 	Uppercase ID for the dist function to use
# returns 	The distance as a % between sequences in an alignFile
def computeDist(alignFile, distFn):
	distFns = {"OUTTER_GAP_CONSERVED":computeDist_outterGapConserved}
	if distFn in distFns.keys():
		return distFns[distFn](alignFile)
	else:
		raise NameError("Error: Invalid grading Fn.  Valid options are:" + distFns.keys())

# t {{}{}} 		A tree from loadTree()
# s {{}{}} 		A tree form loadSequences()
# tax a;b;c;d 		A semicolon-delimited taxonomy string
# upToLvl #		A digit between 1 and 8 indicating depth of tax
# refDB "" 		File path to the fasta DB
# hitsFastaFilePath ""	Output name for FASTA of matching Taxonomic entries and their sequences
# distFn ""		Uppercase ID for the dist function to use
def computeConservPct(t, s, tax, upToLvl, refDB, hitsFastaFilePath, distFn):
	child_seqs = getSequencesOf(t, s, tax, upToLvl)
	print list(chain.from_iterable(child_seqs))
	print("Searching Fasta...")
	children = searchFastaByID(refDB, chain.from_iterable(child_seqs))
	SeqIO.write(children, hitsFastaFilePath, "fasta")
	alignFile = clustalAlign(hitsFastaFilePath)
	pct = computeDist(alignFile, distFn)
	return pct


#Load Tree and Sequences
'''
t=loadTree("data/Unique_taxonomy.lines")
s=loadSequences("data/seq_lin.mapping")
DB="data/bold_unique.fasta"
outFile="temp.fasta"
tax="Animalia;Arthropoda;Remipedia;Nectiopoda"
upToLvl=8
distFn="OUTTER_GAP_CONSERVED"
computeConservPct(t, s, tax, upToLvl, DB, outFile, distFn)
'''
hitsFastaFilePath="data/test.fasta"
distFn="OUTTER_GAP_CONSERVED"
alignFile = "data/test.fasta.aln"#clustalAlign(hitsFastaFilePath)
pct = computeDist(alignFile, distFn)
print(pct)


