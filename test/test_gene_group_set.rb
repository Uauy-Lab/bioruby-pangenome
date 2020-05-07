require "test/unit"
gempath = File.dirname(File.dirname(__FILE__))
path = gempath + '/lib/bio-pangenome.rb'
require path

 
class Test::PangenomeGeneGroup < Test::Unit::TestCase

	def setup
		@genes = Bio::Pangenome.load_genes("./test/data/genes.txt")
		@lines = Bio::Pangenome.load_lines("./test/data/lines.txt")
		@projected_genes = Bio::Pangenome.load_projected_genes("./test/data/mapping.csv.gz", genes: @genes)
		@sequences = Bio::Pangenome.load_sequences(varieties: @lines, genes: @projected_genes, gene_order: @genes, prefix: "./test/data/2000bp", suffix: "_2000bp_RefSeqv1.1.fa.gz", set_id: "2000bp")

	end

	def teardown
		## Nothing really
	end


	def test_complete_alignments
		complete_genes = @sequences.complete
		assert_equal(18, complete_genes.length)
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
		assert(complete_genes.keys.include? "TraesCS6A02G077500")
	end


	def test_complete_haplotype
		complete_genes = @sequences.complete
	
		File.open("./test/data/out/haps.fa", "w") do |out|	
			complete_genes.varieties.each do |v|
				out.puts ">#{v}\n#{complete_genes.haplotype(v,sep: "|")}"
			end
		end

		File.open("./test/data/out/haps.csv", "w") do |out|
			complete_genes.haplotype_matrix.each do |row|
				#puts row.inspect
				out.puts row.join(",")
			end
		end	

		File.open("./test/data/out/haps_all.fa", "w") do |out|	
			complete_genes.varieties.each do |v|
				out.puts ">#{v}\n#{complete_genes.haplotype(v,sep: "|", max_mismatches_per_kbp: 1000)}"
			end
		end

		File.open("./test/data/out/haps_all.csv", "w") do |out|
			complete_genes.haplotype_matrix(max_mismatches_per_kbp: 1000).each do |row|
				#puts row.inspect
				out.puts row.join(",")
			end
		end	

	end
end