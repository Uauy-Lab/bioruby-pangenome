require "test/unit"
gempath = File.dirname(File.dirname(__FILE__))
path = gempath + '/lib/bio-pangenome.rb'
require path

 
class Test::GeneGroupSetHaplotypes < Test::Unit::TestCase

	def setup
		#@genes = Bio::Pangenome.load_genes("./test/data/genes_reg.txt")
		#@lines = Bio::Pangenome.load_lines("./test/data/lines.txt")
		#@projected_genes = Bio::Pangenome.load_projected_genes("./test/data/mapping_reg.csv.gz", genes: @genes)
		#@sequences = Bio::Pangenome.load_sequences(varieties: @lines, genes: @projected_genes, gene_order: @genes, prefix: "./test/data/reg/2000bp", suffix: "_2000bp_RefSeqv1.1.fa.gz", set_id: "2000bp")
		#@complete_genes = @sequences.complete
		
		#puts @sequences.keys
		

	end

	def teardown
		## Nothing really
	end

	def test_load_haplotype_groups
		all_haps = Bio::Pangenome.load_haplotye_groups("./test/data/genes_reg.txt", "./test/data/lines_reg.txt")

		assert_equal(15, all_haps.varieties.size)
		assert(all_haps.blocks.include? "G1" )
		assert(all_haps.blocks.include? "G2" )
		assert_equal(2, all_haps.blocks.size)

		g1 = all_haps.blocks["G1"]
		g2 = all_haps.blocks["G2"]

		assert_equal("CS", g1.variety_haplotype["chinese"])
		assert_equal("P1", g1.variety_haplotype["cadenza"])
		assert_equal("P1", g1.variety_haplotype["paragon"])
		assert_equal("P1", g1.variety_haplotype["norin61"])
		assert_equal("P1", g1.variety_haplotype["lancer"])
		assert_equal("P1", g1.variety_haplotype["claire"])
		assert_equal("P1", g1.variety_haplotype["jagger"])
		assert_equal("P1", g1.variety_haplotype["sy_mattis"])
		assert_equal("P2", g1.variety_haplotype["robigus"])
		assert_equal("P2", g1.variety_haplotype["arinalrfor"])
		assert_equal("P3", g1.variety_haplotype["julius"])
		assert_equal("P4", g1.variety_haplotype["weebil"])
		assert_equal("P5", g1.variety_haplotype["landmark"])
		assert_equal("P6", g1.variety_haplotype["mace"])
		assert_equal("P6", g1.variety_haplotype["stanley"])

		assert(g1.varieties_for["CS"].include?  "chinese")
		assert_equal(g1.varieties_for["P1"], ["cadenza", "paragon", "norin61", "lancer","claire", "jagger", "sy_mattis"])
		assert_equal(g1.varieties_for["P2"], ["robigus","arinalrfor"])
		assert_equal(g1.varieties_for["P3"], ["julius"])
		assert_equal(g1.varieties_for["P4"], ["weebil"])
		assert_equal(g1.varieties_for["P5"], ["landmark"])
		assert_equal(g1.varieties_for["P6"], ["mace", "stanley"])
		assert_raise(FrozenError){ g1.variety_haplotype["chinese"] = "P7" }
	 	assert_raise(FrozenError){ g1.variety_haplotype["P7"]      = "P7" }

	 	assert_equal(["TraesCS6A02G182400",
				"TraesCS6A02G182700",
				"TraesCS6A02G185600",
				"TraesCS6A02G189600",
				"TraesCS6A02G194200",
				"TraesCS6A02G194600",
				"TraesCS6A02G194700",
				"TraesCS6A02G196000",
				"TraesCS6A02G197200",
				"TraesCS6A02G198400"], g1.genes)

		assert_equal(["TraesCS6A02G199300",
			 	"TraesCS6A02G199800",
				"TraesCS6A02G201400",
				"TraesCS6A02G204800",
				"TraesCS6A02G209000",
				"TraesCS6A02G209200",
				"TraesCS6A02G210700",
				"TraesCS6A02G210900",
				"TraesCS6A02G212300",
				"TraesCS6A02G213500",
				"TraesCS6A02G214600",
				"TraesCS6A02G217200"], g2.genes)
	end


	def test_load_haplotypes_groups_chunks
		b1 = Bio::Pangenome.load_haplotye_groups("./test/data/genes_reg.txt", "./test/data/lines_reg.txt", chunk_no: 0, chunks: 2)
		b2 = Bio::Pangenome.load_haplotye_groups("./test/data/genes_reg.txt", "./test/data/lines_reg.txt", chunk_no: 1, chunks: 2) 

		assert_equal(["TraesCS6A02G182400",
				"TraesCS6A02G182700",
				"TraesCS6A02G185600",
				"TraesCS6A02G189600",
				"TraesCS6A02G194200",
				"TraesCS6A02G194600",
				"TraesCS6A02G194700",
				"TraesCS6A02G196000",
				"TraesCS6A02G197200",
				"TraesCS6A02G198400"], b1.genes)

		assert_equal(["TraesCS6A02G199300",
			 	"TraesCS6A02G199800",
				"TraesCS6A02G201400",
				"TraesCS6A02G204800",
				"TraesCS6A02G209000",
				"TraesCS6A02G209200",
				"TraesCS6A02G210700",
				"TraesCS6A02G210900",
				"TraesCS6A02G212300",
				"TraesCS6A02G213500",
				"TraesCS6A02G214600",
				"TraesCS6A02G217200"], b2.genes)

	end


	# def test_aligned_sequences	
	# 	stats_out = File.open("./test/data/out/aligned_stats_reg.csv", "w")	
	# 	stats_keys = [:gaps, :mismatches, :matches, :length, :mismatches_per_kbp, :missing]
	# 	stats_out.puts "gene," + stats_keys.map{|k| k.to_s}.join(",")
	# 	@sequences.complete.each_pair do |name, val|
	# 		seqs = val.aligned_sequences
	# 		mask = val.mask
	# 		File.open("./test/data/out/aligned_#{name}.fa", "w") do |out|	
	# 			out.write mask.output_fasta
	# 		end
	# 		mm = ""
	# 		mm += mask.snp_positions.join(":") if mask.stats[:mismatches_per_kbp] < 5
	# 		stats_out.puts  name + "," + stats_keys.map{|k| mask.stats[k]}.join(",") + "," + mm 
	# 	end
	# 	stats_out.close
	# end

	# def test_complete_haplotype
	# 	complete_genes = @sequences.complete
	# 	#puts complete_genes.inspect
	# 	File.open("./test/data/out/haps_reg.fa", "w") do |out|	
	# 		complete_genes.varieties.each do |v|
	# 			out.puts ">#{v}\n#{complete_genes.haplotype(v,sep: "|")}"
	# 		end
	# 	end

	# 	File.open("./test/data/out/haps_reg.csv", "w") do |out|
	# 		complete_genes.haplotype_matrix.each do |row|
	# 			puts row.inspect
	# 			out.puts row.join(",")
	# 		end
	# 	end	
	# end
end


