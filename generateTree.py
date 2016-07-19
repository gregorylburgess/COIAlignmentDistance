import collections
from Bio import SeqIO


def Tree():
    return collections.defaultdict(Tree)

def loadTree(HeaderFileName):
    t = Tree()
    for line in open(HeaderFileName):
        annotation = line.rstrip()
        params = [x if x else "Unknown" for x in annotation.split(";")]
        add(t, params)
    return t

def loadSequences(SequencesFileName):
    # return a dictionary where each lineage is associated with a list of sequences
    lineageToSequences = collections.defaultdict(list)
    for line in open(SequencesFileName):
        line = line.rstrip()
        data = line.split("\t")
        lineageToSequences[data[1]].append(data[0])
    return lineageToSequences

def add(t, path):
  for node in path:
    t = t[node]

def getSubTree (t, tax):
    myDict = t[tax.split(";")[0]]
    if tax[-1] == ";":
        tax = tax[:-1]
    for taxId in tax.split(";")[1:]:
        myDict =  myDict[taxId]
    return myDict


lineages = []

def getLineagesOf(t, tax, upToLvl=8):
	'''Lists the children OTU of tax, of depth/length upToLvl.

		t	{} Tree Dictionary
		tax	"" Taxonomic name delimited by semicolons.	
		upToLvl # the classification level of the returned items.

		return	List of OTUs of depth upToLvl.
	'''
	global sequenceIds
	lineages=[]
	def getChildren(t, tax, upToLvl=8):
		global sequenceIds
		if tax[-1] == ";":
		    tax = tax[:-1]
		for child in getSubTree(t, tax).keys():
		    getChildren(t, tax+";"+child, upToLvl)

		    if len((tax+";"+child).split(";")) ==  upToLvl:
		        #print tax+";"+child, len((tax+";"+child).split(";"))
		        lineages.append((tax+";"+child).replace("Unknown", ""))
	getChildren(t, tax, upToLvl)
	return lineages

sequenceIds = []
## getChildrenOf(t, lineageToSequences, "Animalia;Arthropoda;Remipedia;Nectiopoda", 8)
def getSequencesOf(t, lineageToSequences, tax, upToLvl=8):
    global sequenceIds
    sequenceIds=[]
    def getChildren(t, lineageToSequences, tax, upToLvl=8):
        global sequenceIds
        if tax[-1] == ";":
            tax = tax[:-1]
        for child in getSubTree(t, tax).keys():
            getChildren(t, lineageToSequences, tax+";"+child, upToLvl)

            if len((tax+";"+child).split(";")) ==  upToLvl:
                #print tax+";"+child, len((tax+";"+child).split(";"))
                sequenceIds.append(lineageToSequences[(tax+";"+child).replace("Unknown", "")])

    getChildren(t, lineageToSequences, tax, upToLvl)
    return sequenceIds


        
