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
		#puts @sequences['TraesCS6A02G081800'].inspect
		assert_equal(2, @sequences['TraesCS6A02G002800'].lengths.uniq.length)
		assert_equal(4, @sequences['TraesCS6A02G006600'].lengths.uniq.length)
		assert_equal(3, @sequences['TraesCS6A02G008400'].lengths.uniq.length)
		assert_equal(2, @sequences['TraesCS6A02G017400'].lengths.uniq.length)
		assert_equal(2, @sequences['TraesCS6A02G027000'].lengths.uniq.length)
		assert_equal(4, @sequences['TraesCS6A02G032300'].lengths.uniq.length)
		assert_equal(4, @sequences['TraesCS6A02G052900'].lengths.uniq.length)
		assert_equal(1, @sequences['TraesCS6A02G077500'].lengths.uniq.length)
		assert_equal(1, @sequences['TraesCS6A02G081800'].lengths.uniq.length)
		assert_equal(14, @sequences['TraesCS6A02G084500'].lengths.uniq.length)
		assert_equal(1, @sequences['TraesCS6A02G096900'].lengths.uniq.length)
		assert_equal(1, @sequences['TraesCS6A02G097400'].lengths.uniq.length)
		assert_equal(1, @sequences['TraesCS6A02G103000'].lengths.uniq.length)
		assert_equal(1, @sequences['TraesCS6A02G120500'].lengths.uniq.length)
		assert_equal(1, @sequences['TraesCS6A02G134500'].lengths.uniq.length)
		assert_equal(3, @sequences['TraesCS6A02G137400'].lengths.uniq.length)
		assert_equal(2, @sequences['TraesCS6A02G143300'].lengths.uniq.length)
		assert_equal(1, @sequences['TraesCS6A02G163100'].lengths.uniq.length)
		assert_equal(4, @sequences['TraesCS6A02G177100'].lengths.uniq.length)
		assert_equal(3, @sequences['TraesCS6A02G177700'].lengths.uniq.length)
		assert_equal(14, @sequences['TraesCS6A02G183800'].lengths.uniq.length)
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
		assert_equal(2, @sequences['TraesCS6A02G307900'].lengths.uniq.length)
		assert_equal(1, @sequences['TraesCS6A02G310500'].lengths.uniq.length)
		assert_equal(2, @sequences['TraesCS6A02G312900'].lengths.uniq.length)
		assert_equal(13, @sequences['TraesCS6A02G329200'].lengths.uniq.length)
		assert_equal(2, @sequences['TraesCS6A02G329600'].lengths.uniq.length)
		assert_equal(1, @sequences['TraesCS6A02G339000'].lengths.uniq.length)
		assert_equal(3, @sequences['TraesCS6A02G389900'].lengths.uniq.length)
		assert_equal(4, @sequences['TraesCS6A02G399100'].lengths.uniq.length)
		assert_equal(3, @sequences['TraesCS6A02G403300'].lengths.uniq.length)
		assert_equal(2, @sequences['TraesCS6A02G012900'].lengths.uniq.length)
		assert_equal(8, @sequences['TraesCS6A02G087800'].lengths.uniq.length)
		assert_equal(6, @sequences['TraesCS6A02G020700'].lengths.uniq.length)
		assert_equal(3, @sequences['TraesCS6A02G025000'].lengths.uniq.length)
		assert_equal(2, @sequences['TraesCS6A02G090800'].lengths.uniq.length)
		assert_equal(1, @sequences['TraesCS6A02G206300'].lengths.uniq.length)
		assert_equal(2, @sequences['TraesCS6A02G378300'].lengths.uniq.length)

		@sequences.each_pair do |name, val|  
			File.open("./test/data/out/#{name}.fa", "w") do |out|
				seqs = val.sequences
				out.write seqs.output_fasta
			end
		end

	end

	def test_aligned_sequences	
		stats_out = File.open("./test/data/out/aligned_stats.csv", "w")	
		stats_keys = [:gaps, :mismatches, :matches, :length, :mismatches_per_kbp, :missing]
		stats_out.puts "gene," + stats_keys.map{|k| k.to_s}.join(",")
		@sequences.complete.each_pair do |name, val|
			seqs = val.aligned_sequences
			mask = val.mask
			File.open("./test/data/out/aligned_#{name}.fa", "w") do |out|	
				out.write mask.output_fasta
			end
			mm = ""
			mm += mask.snp_positions.join(":") if mask.stats[:mismatches_per_kbp] < 5
			stats_out.puts  name + "," + stats_keys.map{|k| mask.stats[k]}.join(",") + "," + mm 
		end
		stats_out.close

		assert_equal(1,   @sequences["TraesCS6A02G256700"].mask.valid_snps.size)
		assert_equal(735, @sequences["TraesCS6A02G256700"].mask.valid_snps[0] )
		positions = Hash.new
		positions['TraesCS6A02G027000'] = [663,1455,1522,3445,3596,3954,4056,4384,4495,4568]
		positions['TraesCS6A02G077500'] = [2821,4030,4403,4773]
		positions['TraesCS6A02G084500'] = [63,117,382,995,1720,1932,2648,3105,4909,4965,5155,5360,5659,5812,6359,6478,6635,8010,8204,8383]
		positions['TraesCS6A02G134500'] = [263,330,986,1238,1793,2103,2488,3393,4268,4360]
		positions['TraesCS6A02G137400'] = [357,1185,1703,2695,2878,3291,3360]
		positions['TraesCS6A02G143300'] = []
		positions['TraesCS6A02G177700'] = [650,1583,5525,6074]
		positions['TraesCS6A02G212100'] = [2837,3091,4679,10090]
		positions['TraesCS6A02G224000'] = [1117,1292,1555,2010,2135,4217,4418,4561]
		positions['TraesCS6A02G235700'] = [193,486,563,950,1997,2214,2525,3089,3335,4031,4293,4910,5054,5111]
		positions['TraesCS6A02G239100'] = [460,1124,2479,3053,4870,5324,5888,6427,6747,7829,8282,8465]
		positions['TraesCS6A02G256700'] = [735]
		positions['TraesCS6A02G261200'] = [645,1926,2136,2752,2852,3046,4084,4207,5848,5986,6642,7048,7153,7267,7853,7949,8043,8121,8183,8256]
		positions['TraesCS6A02G295200'] = [184,469,548,910,1388,1865,3687,3990,4178]
		positions['TraesCS6A02G307900'] = [424,565,722,1145,1387,2004,2321,3464,3598,4434]
		positions['TraesCS6A02G329600'] = []
		positions['TraesCS6A02G389900'] = [993,1567,2153,2290,2742,2856,2913,3190,3429]
		positions['TraesCS6A02G399100'] = [400,1086,1239,1619,1760,1876,1944,3452,3600,3763,3869,4608,5224,5809,5895,6208,6338,6438,6713,6818,7626]

		positions.each_pair do |k,v|
			valid_snps =  @sequences[k].mask.valid_snps
			assert(valid_snps.to_set == v.to_set, "#{k}: #{valid_snps} != #{v}"  )
		end
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
	end

	def test_contigous_blocks
		gene = @sequences["TraesCS6A02G239100"]
	end
end