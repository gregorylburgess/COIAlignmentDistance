from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, assert_true
from generateTree import *

# Simple test taxa to use throughout testing
def simpleTaxPath():
	testTax = ["a;b;c", "x;y;z", "a;d;e"]
	paths =[]
	for tax in testTax:
		paths.append(tax.split(";"))
	return paths

class TestGenerateTree(object):
	# Structure for a simple test tree

	# 	root
	#	a		x	
	#	b	d	y
	#	c	e	z
	@classmethod
	def setup_class(klass):
		"""This method is run once for each class before any tests are run"""

	@classmethod
	def teardown_class(klass):
		"""This method is run once for each class _after_ all tests are run"""
	
	def setUp(self):
		"""This method is run once before _each_ test method is executed"""
		

	def teardown(self):
		"""This method is run once after _each_ test method is executed"""

	# Check that all lines in a file are loaded
	#	 loadTree(HeaderFileName)
	def test_loadTree(self):
		testFileDir="tests/"
		testFiles = ["testData/1lines.lines",
			     "testData/2lines.lines",
			     "testData/5lines.lines",
			     "testData/emptyFile"
			     #"data/Unique_taxonomy.lines"
			    ]
		for testFile in testFiles:
			testFile = testFileDir + testFile
			t=loadTree(testFile)
			lineCount = sum((1 for i in open(testFile, 'rb')))
			print(len(t))
			print(lineCount)
			# Make sure we have didn't miss any lines
			assert_equal(len(t), lineCount)

	# Check that we can get the right subtree.
	#	 getSubTree(t, tax)
	def test_getSubTree(self):
		tree=Tree()
		paths = simpleTaxPath()
		for path in paths:
			add(tree, path)
		for path in paths:
			root=tree
			for node in path:
				assert_true(node in root.keys())
				root=root[node]
		
	# Check that added nodes go to the right parents.
	#	 add(t, path)
	def test_add(self):
		tree=Tree()
		paths = simpleTaxPath()
		paths = simpleTaxPath()
		for path in paths:
			add(tree, path)
		for path in paths:
			root=tree
			for node in path:
				assert_true(node in root.keys())
				root=root[node]
	

	def test_loadSequences(self):
		'''
		s=loadSequences("data/seq_lin.mapping")
		loadSequences(SequencesFileName):
	
	
		fastaFilePath="data/bold_coi_11_05_2015.fasta"
		outFile="temp.fasta"
		tax="Animalia;Arthropoda;Remipedia;Nectiopoda"
		upToLvl=8
		computeConservPct(t, s, tax, upToLvl, fastaFilePath, outFile)

		'''

	'''

	def test_getLineagesOf(self):
		#getLineagesOf(t, tax, upToLvl=8):


	def test_getSequencesOf(self):
		#getSequencesOf(t, lineageToSequences, tax, upToLvl=8):



	def test_run(self):
		a = A()
		assert_equal(a.value, "Some Value")
		assert_not_equal(a.value, "Incorrect Value")

	def test_return_true(self):
		a = A()
		assert_equal(a.return_true(), True)
		assert_not_equal(a.return_true(), False)

	def test_raise_exc(self):
		a = A()
		assert_raises(KeyError, a.raise_exc, "A value")

	@raises(KeyError)
	def test_raise_exc_with_decorator(self):
		a = A()
		a.raise_exc("A message")
	'''
