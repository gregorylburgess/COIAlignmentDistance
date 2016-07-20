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
        data[1] = ";".join([x if x else "Unknown" for x in data[1].split(";")])
        lineageToSequences[data[1]].append(data[0])
    return lineageToSequences

def add(t, path):
  for node in path:
    t = t[node]

def getSubTree (t, tax):
    myDict = t[tax.split(";")[0]]
    while tax[-1] == ";":
        tax = tax[:-1]
    for taxId in tax.split(";")[1:]:
        myDict =  myDict[taxId]
    return myDict


lineages = []

def getLineagesOf(t, tax, upToLvl=8):
    if ";;" in tax:
        print "I am here----"
        raise Exception("empty levels, such as ;; are not permitted in a taxonomy")

    global sequenceIds
    lineages=[]
    def getChildren(t, tax, upToLvl=8):
        global sequenceIds
        while tax[-1] == ";":
            tax = tax[:-1]
        for child in getSubTree(t, tax).keys():
            getChildren(t, tax+";"+child, upToLvl)
            if len((tax+";"+child).split(";")) ==  upToLvl:
                lineages.append((tax+";"+child))
    getChildren(t, tax, upToLvl)
    return lineages

sequenceIds = []
def getSequencesOf(t, lineageToSequences, tax, upToLvl=8):
    global sequenceIds
    sequenceIds=[]
    def getChildren(t, lineageToSequences, tax, upToLvl=8):
        while tax[-1] == ";":
            tax = tax[:-1]
        if len(tax.split(";")) ==  upToLvl:
            #print "Getting seqeuence of %s" % tax
            sequenceIds.append(lineageToSequences[tax])
        else:
            for child in getSubTree(t, tax).keys():
                getChildren(t, lineageToSequences, tax+";"+child, upToLvl)

    getChildren(t, lineageToSequences, tax, upToLvl)
    return sequenceIds


        

