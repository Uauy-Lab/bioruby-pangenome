require "damerau-levenshtein"
module Bio::Pangenome
	class GeneGroupSet < Hash 
		attr_accessor :varieties
		attr_accessor :genes 
		attr_accessor :variety_haplotype
		def complete
			ret = GeneGroupSet.new { |h, k| h[k] = GeneGroup.new }
			ret.varieties = varieties
			ret.genes = []
			self.each_pair do |k,v|
				ret[k] = v if v.keys.length == varieties.length
			end 
			ret.genes = genes.select {|g| ret.include? g}
			return ret
		end

		def varieties_for
			return @varieties_for if @varieties_for
			@variety_haplotype.freeze
			@varieties_for = Hash.new{|h,k| h[k] = []}
			variety_haplotype.each_pair do |name, val| 
				@varieties_for[val] << name 
			end
			return @varieties_for
		end

		def haplotype(variety, sep: "", max_mismatches_per_kbp: 5)
			hap = genes.map do |gene|
				next if self[gene].mask.stats[:mismatches_per_kbp] > max_mismatches_per_kbp
				self[gene].snps_hap[variety]
			end
			hap.select{|h| not h.nil?  and h.size > 0}.join(sep)
		end

		def haplotypes(max_mismatches_per_kbp: 5)
			varieties.map{|v| haplotype(v) }
		end
 
		def haplotype_groups(max_mismatches_per_kbp: 5)
			hap_mat = haplotypes(max_mismatches_per_kbp: max_mismatches_per_kbp)
			relations = []
			varieties.each_with_index do |v1,i|
				varieties.slice(0, i).each_with_index do |v2, j|
					dist = DamerauLevenshtein.distance(hap_mat[i], hap_mat[j])
					relations << Bio::Relation.new(v1, v2, dist)
				end
			end
			dist_mat = Bio::Pathway.new(relations, 'undirected')
			puts dist_mat
			dist_mat
		end

		def haplotype_matrix(max_mismatches_per_kbp: 5)
			ret = []
			tmp = ["gene", "position"]
			tmp.concat varieties
			ret << tmp
			genes.each do |gene|
				next if self[gene].mask.stats[:mismatches_per_kbp] > max_mismatches_per_kbp
				ret.concat self[gene].haplotype_matrix
			end
			ret
		end
	end

	def self.load_sequences( varieties:[], genes:{}, gene_order:[], prefix: "../flanking/filtered/",  suffix: ".cds.fa.gz", set_id: "cds" )
		ret = GeneGroupSet.new { |h, k| h[k] = GeneGroup.new }
		ret.varieties = varieties
		ret.genes = gene_order
		varieties.each do |variety|
			path = "#{prefix}/#{variety}#{suffix}"
			Bio::Extensions.fasta_gz_iterator(path) do |entry|
				definition, region = entry.definition.split("::")
				seq_name = parseSequenceName(region, definition) if region
				seq_name = GeneFlankingRegion.new(entry.definition, nil, "",
					"", entry.definition, set_id, nil, variety ) unless region
				arr = definition.split(".")
				next unless genes.include? seq_name.gene
				row = genes[seq_name.gene]
				row = genes[arr[0]] unless region
				seq_name.gene = row["gene"]
				seq = entry.seq
				seq_name.sequence = seq
				base_gene = row["gene"]
				ret[base_gene][variety] = seq_name unless ret[base_gene][variety]
				ret[base_gene].gene = base_gene
				ret[base_gene].expected_varieties = varieties
			end
		end
		ret
	end

	class GeneGroupSetHaplotypes 
		attr_reader :varieties, :all_genes, :chunks, :chunk_no, :blocks
		def initialize(chunks: 1, chunk_no: 0)
			@varieties = Set.new 
			@all_genes = Set.new
			@chunks = chunks
			@chunk_no = chunk_no
			@blocks = Hash.new do|h,k| 
				tmp = GeneGroupSet.new 
				tmp.variety_haplotype = Hash.new  
				tmp.genes     = []
				tmp.varieties = Set.new
				h[k]=tmp
			end
		end
	end

	def self.load_haplotye_groups(genes_file, haplotypes_file, chunk_no: 0, chunks: 1)
		ret = GeneGroupSetHaplotypes.new(chunk_no: chunk_no, chunks: chunks)
		CSV.foreach(haplotypes_file, col_sep:"\t", headers: true) do |row|
			ret.varieties << row["line"]
			ret.blocks[row["block"]].varieties << row["line"]
			ret.blocks[row["block"]].variety_haplotype[row["line"]] = row["haplotype"]
		end


		
		CSV.foreach(genes_file, col_sep:"\t", headers: true) do |row|
			ret.all_genes << row["gene"]
			ret.blocks[row["block"]].genes << row["gene"]
		end
		#puts ret.inspect
		ret
	end
end