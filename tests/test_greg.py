from nose.tools import assert_equal
from nose.tools import assert_not_equal
from nose.tools import assert_raises
from nose.tools import raises
import greg

class TestA(object):

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

	def test_parseFASTA(self):
		return True
		#parse empty file
		#parse single entry file
		#parse multiple entry file

	def test_clustalAlign(self):
		return True
		#align empty file
		#align single entry file
		#align multiple entry file

	def test_parseALN(self):
		return True
		#parse empty file
		#parse single entry file
		#parse multiple entry file

	def test_searchFastabyID(self):
		return True
		#search each k,p,c,o,f,g,s for existing sp
		#search each k,p,c,o,f,g,s for non-existant

	def test_gradeAlignment_outterGapConserved(self):
		#test 
		return True
	def test_defGradeAlignmnet(self):
		#test existing alg
		#test non-existing alg
		#test 		
		return True
	def computeConservPct(self):
		# Load Tree and Sequences
		'''
		t=loadTree("data/Unique_taxonomy.lines")
		s=loadSequences("data/seq_lin.mapping")
		fastaFilePath="data/bold_coi_11_05_2015.fasta"
		outFile="temp.fasta"
		tax="Animalia;Arthropoda;Remipedia;Nectiopoda"
		upToLvl=8
		computeConservPct(t, s, tax, upToLvl, fastaFilePath, outFile)

		'''





