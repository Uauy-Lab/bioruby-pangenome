require "test/unit"
gempath = File.dirname(File.dirname(__FILE__))
path = gempath + '/lib/bio-pangenome.rb'
require path

 
class Test::PangenomeGeneGroup < Test::Unit::TestCase
	def setup
		@genes = Bio::Pangenome.load_genes("./test/data/genes.txt")
		@lines = Bio::Pangenome.load_lines("./test/data/lines.txt")
		@projected_genes = Bio::Pangenome.load_projected_genes("./test/data/mapping.csv.gz", genes: @genes)
		@sequences = Bio::Pangenome.load_sequences(varieties: @lines, genes: @projected_genes, prefix: "./test/data/2000bp", suffix: "_2000bp_RefSeqv1.1.fa.gz", set_id: "2000bp")

		@gene_1 =  @sequences["TraesCS6A02G329200"]
		@gene_2 =  @sequences["TraesCS6A02G310500"]

	end

	def teardown
		## Nothing really
	end

	def test_sequences
		assert_equal(2, @sequences['TraesCS6A02G002800'].lengths.uniq.length)
		assert_equal(4, @sequences['TraesCS6A02G006600'].lengths.uniq.length)
		assert_equal(3, @sequences['TraesCS6A02G008400'].lengths.uniq.length)
		assert_equal(2, @sequences['TraesCS6A02G017400'].lengths.uniq.length)
		assert_equal(2, @sequences['TraesCS6A02G027000'].lengths.uniq.length)
		assert_equal(4, @sequences['TraesCS6A02G032300'].lengths.uniq.length)
		assert_equal(4, @sequences['TraesCS6A02G052900'].lengths.uniq.length)
		assert_equal(1, @sequences['TraesCS6A02G077500'].lengths.uniq.length)
		assert_equal(6, @sequences['TraesCS6A02G081800'].lengths.uniq.length)
		assert_equal(14, @sequences['TraesCS6A02G084500'].lengths.uniq.length)
		assert_equal(2, @sequences['TraesCS6A02G096900'].lengths.uniq.length)
		assert_equal(1, @sequences['TraesCS6A02G097400'].lengths.uniq.length)
		assert_equal(1, @sequences['TraesCS6A02G103000'].lengths.uniq.length)
		assert_equal(1, @sequences['TraesCS6A02G120500'].lengths.uniq.length)
		assert_equal(1, @sequences['TraesCS6A02G134500'].lengths.uniq.length)
		assert_equal(5, @sequences['TraesCS6A02G137400'].lengths.uniq.length)
		assert_equal(3, @sequences['TraesCS6A02G143300'].lengths.uniq.length)
		assert_equal(1, @sequences['TraesCS6A02G163100'].lengths.uniq.length)
		assert_equal(4, @sequences['TraesCS6A02G177100'].lengths.uniq.length)
		assert_equal(3, @sequences['TraesCS6A02G177700'].lengths.uniq.length)
		assert_equal(15, @sequences['TraesCS6A02G183800'].lengths.uniq.length)
		assert_equal(2, @sequences['TraesCS6A02G197000'].lengths.uniq.length)
		assert_equal(2, @sequences['TraesCS6A02G212100'].lengths.uniq.length)
		assert_equal(1, @sequences['TraesCS6A02G224000'].lengths.uniq.length)
		assert_equal(2, @sequences['TraesCS6A02G235700'].lengths.uniq.length)
		assert_equal(1, @sequences['TraesCS6A02G239100'].lengths.uniq.length)
		assert_equal(1, @sequences['TraesCS6A02G243900'].lengths.uniq.length)
		assert_equal(1, @sequences['TraesCS6A02G256700'].lengths.uniq.length)
		assert_equal(2, @sequences['TraesCS6A02G261200'].lengths.uniq.length)
		assert_equal(1, @sequences['TraesCS6A02G264300'].lengths.uniq.length)
		assert_equal(4, @sequences['TraesCS6A02G418100'].lengths.uniq.length)
		assert_equal(1, @sequences['TraesCS6A02G282100'].lengths.uniq.length)
		assert_equal(2, @sequences['TraesCS6A02G295200'].lengths.uniq.length)
		assert_equal(4, @sequences['TraesCS6A02G307900'].lengths.uniq.length)
		assert_equal(1, @sequences['TraesCS6A02G310500'].lengths.uniq.length)
		assert_equal(7, @sequences['TraesCS6A02G312900'].lengths.uniq.length)
		assert_equal(14, @sequences['TraesCS6A02G329200'].lengths.uniq.length)
		assert_equal(2, @sequences['TraesCS6A02G329600'].lengths.uniq.length)
		assert_equal(1, @sequences['TraesCS6A02G339000'].lengths.uniq.length)
		assert_equal(3, @sequences['TraesCS6A02G389900'].lengths.uniq.length)
		assert_equal(4, @sequences['TraesCS6A02G399100'].lengths.uniq.length)
		assert_equal(6, @sequences['TraesCS6A02G403300'].lengths.uniq.length)
		assert_equal(2, @sequences['TraesCS6A02G012900'].lengths.uniq.length)
		assert_equal(8, @sequences['TraesCS6A02G087800'].lengths.uniq.length)
		assert_equal(7, @sequences['TraesCS6A02G020700'].lengths.uniq.length)
		assert_equal(3, @sequences['TraesCS6A02G025000'].lengths.uniq.length)
		assert_equal(2, @sequences['TraesCS6A02G090800'].lengths.uniq.length)
		assert_equal(1, @sequences['TraesCS6A02G206300'].lengths.uniq.length)
		assert_equal(2, @sequences['TraesCS6A02G378300'].lengths.uniq.length)

		@sequences.each_pair do |name, val|  
			File.open("./test/data/out/#{name}.fa", "w") do |out|
				seqs = val.sequences
				out.write seqs.to_fasta
			end
		end

	end

	def test_aligned_sequences		
		@sequences.complete.each_pair do |name, val|  
			File.open("./test/data/out/aligned_#{name}.fa", "w") do |out|
				seqs = val.aligned_sequences
				out.write seqs.to_fasta
			end
		end
	end

	def test_complete_alignments
		complete_genes = @sequences.complete
		assert_equal(16, complete_genes.length)
		assert(complete_genes.keys.include? "TraesCS6A02G027000") 
		assert(complete_genes.keys.include? "TraesCS6A02G134500") 
		assert(complete_genes.keys.include? "TraesCS6A02G137400") 
		assert(complete_genes.keys.include? "TraesCS6A02G143300") 
		assert(complete_genes.keys.include? "TraesCS6A02G177700") 
		assert(complete_genes.keys.include? "TraesCS6A02G212100") 
		assert(complete_genes.keys.include? "TraesCS6A02G224000") 
		assert(complete_genes.keys.include? "TraesCS6A02G235700") 
		assert(complete_genes.keys.include? "TraesCS6A02G239100") 
		assert(complete_genes.keys.include? "TraesCS6A02G256700") 
		assert(complete_genes.keys.include? "TraesCS6A02G261200") 
		assert(complete_genes.keys.include? "TraesCS6A02G295200") 
		assert(complete_genes.keys.include? "TraesCS6A02G307900") 
		assert(complete_genes.keys.include? "TraesCS6A02G329600") 
		assert(complete_genes.keys.include? "TraesCS6A02G389900") 
		assert(complete_genes.keys.include? "TraesCS6A02G399100")
	end

	def test_contigous_blocks
		gene = @sequences["TraesCS6A02G239100"]
	end
end