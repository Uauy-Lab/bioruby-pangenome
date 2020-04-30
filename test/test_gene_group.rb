require "test/unit"
gempath = File.dirname(File.dirname(__FILE__))
path = gempath + '/lib/bio-pangenome.rb'
require path

 
class TestBioPangenomeGeneGroup < Test::Unit::TestCase
	def setup
		@genes = BioPangenome.load_genes("./test/data/genes.txt")
		@lines = BioPangenome.load_lines("./test/data/lines.txt")
		@projected_genes = BioPangenome.load_projected_genes("./test/data/mapping.csv.gz", genes: @genes)
		@sequences = BioPangenome.load_sequences(varieties: @lines, genes: @projected_genes, prefix: "./test/data/2000bp", suffix: "_2000bp_RefSeqv1.1.fa.gz", set_id: "2000bp")

		@gene_1 =  @sequences["TraesCS6A02G329200"]
		@gene_2 =  @sequences["TraesCS6A02G310500"]

	end

	def teardown
		## Nothing really
	end

	def test_aligned_sequences		
		#@gene_1.sequences 
		txt = ""
		@sequences.each_pair do |name, val|  
			lengths = val.lengths
			lengths.uniq!
			txt +=  "assert_equal(#{lengths.size}, @sequences['#{name}'].lengths.uniq!)\n"
		end
		puts txt
	end

end