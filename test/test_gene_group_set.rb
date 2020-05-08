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
		#@complete_genes = @sequences.complete
		
		@region_genes = Bio::Pangenome.load_genes("./test/data/genes_block.txt")
		#@projected_genes_ = Bio::Pangenome.load_projected_genes("./test/data/mapping.csv.gz", genes: @genes)
		@sequences_region = Bio::Pangenome.load_sequences(varieties: @lines, genes: @projected_genes, gene_order: @region_genes, prefix: "./test/data/2000bp", suffix: "_2000bp_RefSeqv1.1.fa.gz", set_id: "2000bp")
		
		

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

	def test_haplotype_groups
		hap_groups = @sequences_region.haplotype_groups
		#puts hap_groups.dump_list()
		#puts hap_groups.graph["arinalrfor"]["cadenza"]
		#hap_groups.graph.each do |from, h|
		#	h.each do |to, relation|
		#		puts "assert_equal(#{relation}, hap_groups['#{from}']['#{to}'])"
		#	end
		#end

		assert_equal(11, hap_groups.graph['cadenza']['arinalrfor'])
		assert_equal(11, hap_groups.graph['cadenza']['chinese'])
		assert_equal(11, hap_groups.graph['cadenza']['claire'])
		assert_equal(11, hap_groups.graph['cadenza']['jagger'])
		assert_equal(11, hap_groups.graph['cadenza']['julius'])
		assert_equal(11, hap_groups.graph['cadenza']['lancer'])
		assert_equal(11, hap_groups.graph['cadenza']['landmark'])
		assert_equal(11, hap_groups.graph['cadenza']['mace'])
		assert_equal(11, hap_groups.graph['cadenza']['norin61'])
		assert_equal(0, hap_groups.graph['cadenza']['paragon'])
		assert_equal(11, hap_groups.graph['cadenza']['robigus'])
		assert_equal(11, hap_groups.graph['cadenza']['stanley'])
		assert_equal(11, hap_groups.graph['cadenza']['sy_mattis'])
		assert_equal(5, hap_groups.graph['cadenza']['weebil'])
		assert_equal(11, hap_groups.graph['arinalrfor']['cadenza'])
		assert_equal(6, hap_groups.graph['arinalrfor']['chinese'])
		assert_equal(6, hap_groups.graph['arinalrfor']['claire'])
		assert_equal(6, hap_groups.graph['arinalrfor']['jagger'])
		assert_equal(1, hap_groups.graph['arinalrfor']['julius'])
		assert_equal(6, hap_groups.graph['arinalrfor']['lancer'])
		assert_equal(0, hap_groups.graph['arinalrfor']['landmark'])
		assert_equal(1, hap_groups.graph['arinalrfor']['mace'])
		assert_equal(6, hap_groups.graph['arinalrfor']['norin61'])
		assert_equal(11, hap_groups.graph['arinalrfor']['paragon'])
		assert_equal(0, hap_groups.graph['arinalrfor']['robigus'])
		assert_equal(1, hap_groups.graph['arinalrfor']['stanley'])
		assert_equal(5, hap_groups.graph['arinalrfor']['sy_mattis'])
		assert_equal(11, hap_groups.graph['arinalrfor']['weebil'])
		assert_equal(6, hap_groups.graph['chinese']['arinalrfor'])
		assert_equal(11, hap_groups.graph['chinese']['cadenza'])
		assert_equal(5, hap_groups.graph['chinese']['claire'])
		assert_equal(5, hap_groups.graph['chinese']['jagger'])
		assert_equal(7, hap_groups.graph['chinese']['julius'])
		assert_equal(5, hap_groups.graph['chinese']['lancer'])
		assert_equal(6, hap_groups.graph['chinese']['landmark'])
		assert_equal(7, hap_groups.graph['chinese']['mace'])
		assert_equal(5, hap_groups.graph['chinese']['norin61'])
		assert_equal(11, hap_groups.graph['chinese']['paragon'])
		assert_equal(6, hap_groups.graph['chinese']['robigus'])
		assert_equal(7, hap_groups.graph['chinese']['stanley'])
		assert_equal(4, hap_groups.graph['chinese']['sy_mattis'])
		assert_equal(11, hap_groups.graph['chinese']['weebil'])
		assert_equal(6, hap_groups.graph['claire']['arinalrfor'])
		assert_equal(11, hap_groups.graph['claire']['cadenza'])
		assert_equal(5, hap_groups.graph['claire']['chinese'])
		assert_equal(0, hap_groups.graph['claire']['jagger'])
		assert_equal(7, hap_groups.graph['claire']['julius'])
		assert_equal(2, hap_groups.graph['claire']['lancer'])
		assert_equal(6, hap_groups.graph['claire']['landmark'])
		assert_equal(7, hap_groups.graph['claire']['mace'])
		assert_equal(2, hap_groups.graph['claire']['norin61'])
		assert_equal(11, hap_groups.graph['claire']['paragon'])
		assert_equal(6, hap_groups.graph['claire']['robigus'])
		assert_equal(7, hap_groups.graph['claire']['stanley'])
		assert_equal(1, hap_groups.graph['claire']['sy_mattis'])
		assert_equal(11, hap_groups.graph['claire']['weebil'])
		assert_equal(6, hap_groups.graph['jagger']['arinalrfor'])
		assert_equal(11, hap_groups.graph['jagger']['cadenza'])
		assert_equal(5, hap_groups.graph['jagger']['chinese'])
		assert_equal(0, hap_groups.graph['jagger']['claire'])
		assert_equal(7, hap_groups.graph['jagger']['julius'])
		assert_equal(2, hap_groups.graph['jagger']['lancer'])
		assert_equal(6, hap_groups.graph['jagger']['landmark'])
		assert_equal(7, hap_groups.graph['jagger']['mace'])
		assert_equal(2, hap_groups.graph['jagger']['norin61'])
		assert_equal(11, hap_groups.graph['jagger']['paragon'])
		assert_equal(6, hap_groups.graph['jagger']['robigus'])
		assert_equal(7, hap_groups.graph['jagger']['stanley'])
		assert_equal(1, hap_groups.graph['jagger']['sy_mattis'])
		assert_equal(11, hap_groups.graph['jagger']['weebil'])
		assert_equal(1, hap_groups.graph['julius']['arinalrfor'])
		assert_equal(11, hap_groups.graph['julius']['cadenza'])
		assert_equal(7, hap_groups.graph['julius']['chinese'])
		assert_equal(7, hap_groups.graph['julius']['claire'])
		assert_equal(7, hap_groups.graph['julius']['jagger'])
		assert_equal(7, hap_groups.graph['julius']['lancer'])
		assert_equal(1, hap_groups.graph['julius']['landmark'])
		assert_equal(2, hap_groups.graph['julius']['mace'])
		assert_equal(7, hap_groups.graph['julius']['norin61'])
		assert_equal(11, hap_groups.graph['julius']['paragon'])
		assert_equal(1, hap_groups.graph['julius']['robigus'])
		assert_equal(2, hap_groups.graph['julius']['stanley'])
		assert_equal(6, hap_groups.graph['julius']['sy_mattis'])
		assert_equal(11, hap_groups.graph['julius']['weebil'])
		assert_equal(6, hap_groups.graph['lancer']['arinalrfor'])
		assert_equal(11, hap_groups.graph['lancer']['cadenza'])
		assert_equal(5, hap_groups.graph['lancer']['chinese'])
		assert_equal(2, hap_groups.graph['lancer']['claire'])
		assert_equal(2, hap_groups.graph['lancer']['jagger'])
		assert_equal(7, hap_groups.graph['lancer']['julius'])
		assert_equal(6, hap_groups.graph['lancer']['landmark'])
		assert_equal(7, hap_groups.graph['lancer']['mace'])
		assert_equal(0, hap_groups.graph['lancer']['norin61'])
		assert_equal(11, hap_groups.graph['lancer']['paragon'])
		assert_equal(6, hap_groups.graph['lancer']['robigus'])
		assert_equal(7, hap_groups.graph['lancer']['stanley'])
		assert_equal(1, hap_groups.graph['lancer']['sy_mattis'])
		assert_equal(11, hap_groups.graph['lancer']['weebil'])
		assert_equal(0, hap_groups.graph['landmark']['arinalrfor'])
		assert_equal(11, hap_groups.graph['landmark']['cadenza'])
		assert_equal(6, hap_groups.graph['landmark']['chinese'])
		assert_equal(6, hap_groups.graph['landmark']['claire'])
		assert_equal(6, hap_groups.graph['landmark']['jagger'])
		assert_equal(1, hap_groups.graph['landmark']['julius'])
		assert_equal(6, hap_groups.graph['landmark']['lancer'])
		assert_equal(1, hap_groups.graph['landmark']['mace'])
		assert_equal(6, hap_groups.graph['landmark']['norin61'])
		assert_equal(11, hap_groups.graph['landmark']['paragon'])
		assert_equal(0, hap_groups.graph['landmark']['robigus'])
		assert_equal(1, hap_groups.graph['landmark']['stanley'])
		assert_equal(5, hap_groups.graph['landmark']['sy_mattis'])
		assert_equal(11, hap_groups.graph['landmark']['weebil'])
		assert_equal(1, hap_groups.graph['mace']['arinalrfor'])
		assert_equal(11, hap_groups.graph['mace']['cadenza'])
		assert_equal(7, hap_groups.graph['mace']['chinese'])
		assert_equal(7, hap_groups.graph['mace']['claire'])
		assert_equal(7, hap_groups.graph['mace']['jagger'])
		assert_equal(2, hap_groups.graph['mace']['julius'])
		assert_equal(7, hap_groups.graph['mace']['lancer'])
		assert_equal(1, hap_groups.graph['mace']['landmark'])
		assert_equal(7, hap_groups.graph['mace']['norin61'])
		assert_equal(11, hap_groups.graph['mace']['paragon'])
		assert_equal(1, hap_groups.graph['mace']['robigus'])
		assert_equal(0, hap_groups.graph['mace']['stanley'])
		assert_equal(6, hap_groups.graph['mace']['sy_mattis'])
		assert_equal(11, hap_groups.graph['mace']['weebil'])
		assert_equal(6, hap_groups.graph['norin61']['arinalrfor'])
		assert_equal(11, hap_groups.graph['norin61']['cadenza'])
		assert_equal(5, hap_groups.graph['norin61']['chinese'])
		assert_equal(2, hap_groups.graph['norin61']['claire'])
		assert_equal(2, hap_groups.graph['norin61']['jagger'])
		assert_equal(7, hap_groups.graph['norin61']['julius'])
		assert_equal(0, hap_groups.graph['norin61']['lancer'])
		assert_equal(6, hap_groups.graph['norin61']['landmark'])
		assert_equal(7, hap_groups.graph['norin61']['mace'])
		assert_equal(11, hap_groups.graph['norin61']['paragon'])
		assert_equal(6, hap_groups.graph['norin61']['robigus'])
		assert_equal(7, hap_groups.graph['norin61']['stanley'])
		assert_equal(1, hap_groups.graph['norin61']['sy_mattis'])
		assert_equal(11, hap_groups.graph['norin61']['weebil'])
		assert_equal(11, hap_groups.graph['paragon']['arinalrfor'])
		assert_equal(0, hap_groups.graph['paragon']['cadenza'])
		assert_equal(11, hap_groups.graph['paragon']['chinese'])
		assert_equal(11, hap_groups.graph['paragon']['claire'])
		assert_equal(11, hap_groups.graph['paragon']['jagger'])
		assert_equal(11, hap_groups.graph['paragon']['julius'])
		assert_equal(11, hap_groups.graph['paragon']['lancer'])
		assert_equal(11, hap_groups.graph['paragon']['landmark'])
		assert_equal(11, hap_groups.graph['paragon']['mace'])
		assert_equal(11, hap_groups.graph['paragon']['norin61'])
		assert_equal(11, hap_groups.graph['paragon']['robigus'])
		assert_equal(11, hap_groups.graph['paragon']['stanley'])
		assert_equal(11, hap_groups.graph['paragon']['sy_mattis'])
		assert_equal(5, hap_groups.graph['paragon']['weebil'])
		assert_equal(0, hap_groups.graph['robigus']['arinalrfor'])
		assert_equal(11, hap_groups.graph['robigus']['cadenza'])
		assert_equal(6, hap_groups.graph['robigus']['chinese'])
		assert_equal(6, hap_groups.graph['robigus']['claire'])
		assert_equal(6, hap_groups.graph['robigus']['jagger'])
		assert_equal(1, hap_groups.graph['robigus']['julius'])
		assert_equal(6, hap_groups.graph['robigus']['lancer'])
		assert_equal(0, hap_groups.graph['robigus']['landmark'])
		assert_equal(1, hap_groups.graph['robigus']['mace'])
		assert_equal(6, hap_groups.graph['robigus']['norin61'])
		assert_equal(11, hap_groups.graph['robigus']['paragon'])
		assert_equal(1, hap_groups.graph['robigus']['stanley'])
		assert_equal(5, hap_groups.graph['robigus']['sy_mattis'])
		assert_equal(11, hap_groups.graph['robigus']['weebil'])
		assert_equal(1, hap_groups.graph['stanley']['arinalrfor'])
		assert_equal(11, hap_groups.graph['stanley']['cadenza'])
		assert_equal(7, hap_groups.graph['stanley']['chinese'])
		assert_equal(7, hap_groups.graph['stanley']['claire'])
		assert_equal(7, hap_groups.graph['stanley']['jagger'])
		assert_equal(2, hap_groups.graph['stanley']['julius'])
		assert_equal(7, hap_groups.graph['stanley']['lancer'])
		assert_equal(1, hap_groups.graph['stanley']['landmark'])
		assert_equal(0, hap_groups.graph['stanley']['mace'])
		assert_equal(7, hap_groups.graph['stanley']['norin61'])
		assert_equal(11, hap_groups.graph['stanley']['paragon'])
		assert_equal(1, hap_groups.graph['stanley']['robigus'])
		assert_equal(6, hap_groups.graph['stanley']['sy_mattis'])
		assert_equal(11, hap_groups.graph['stanley']['weebil'])
		assert_equal(5, hap_groups.graph['sy_mattis']['arinalrfor'])
		assert_equal(11, hap_groups.graph['sy_mattis']['cadenza'])
		assert_equal(4, hap_groups.graph['sy_mattis']['chinese'])
		assert_equal(1, hap_groups.graph['sy_mattis']['claire'])
		assert_equal(1, hap_groups.graph['sy_mattis']['jagger'])
		assert_equal(6, hap_groups.graph['sy_mattis']['julius'])
		assert_equal(1, hap_groups.graph['sy_mattis']['lancer'])
		assert_equal(5, hap_groups.graph['sy_mattis']['landmark'])
		assert_equal(6, hap_groups.graph['sy_mattis']['mace'])
		assert_equal(1, hap_groups.graph['sy_mattis']['norin61'])
		assert_equal(11, hap_groups.graph['sy_mattis']['paragon'])
		assert_equal(5, hap_groups.graph['sy_mattis']['robigus'])
		assert_equal(6, hap_groups.graph['sy_mattis']['stanley'])
		assert_equal(11, hap_groups.graph['sy_mattis']['weebil'])
		assert_equal(11, hap_groups.graph['weebil']['arinalrfor'])
		assert_equal(5, hap_groups.graph['weebil']['cadenza'])
		assert_equal(11, hap_groups.graph['weebil']['chinese'])
		assert_equal(11, hap_groups.graph['weebil']['claire'])
		assert_equal(11, hap_groups.graph['weebil']['jagger'])
		assert_equal(11, hap_groups.graph['weebil']['julius'])
		assert_equal(11, hap_groups.graph['weebil']['lancer'])
		assert_equal(11, hap_groups.graph['weebil']['landmark'])
		assert_equal(11, hap_groups.graph['weebil']['mace'])
		assert_equal(11, hap_groups.graph['weebil']['norin61'])
		assert_equal(5, hap_groups.graph['weebil']['paragon'])
		assert_equal(11, hap_groups.graph['weebil']['robigus'])
		assert_equal(11, hap_groups.graph['weebil']['stanley'])
		assert_equal(11, hap_groups.graph['weebil']['sy_mattis'])
		complete_genes = @sequences_region
		File.open("./test/data/out/haps_reg.fa", "w") do |out|	
			complete_genes.varieties.each do |v|
				out.puts ">#{v}\n#{complete_genes.haplotype(v,sep: "|")}"
			end
		end
	end
end