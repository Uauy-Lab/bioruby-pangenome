require "damerau-levenshtein"
module Bio::Pangenome
	class GeneGroupSet < Hash 
		attr_accessor :varieties
		attr_accessor :genes 
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
			#puts hap_mat
			varieties.each_with_index do |v1,i|
				varieties.slice(0, i).each_with_index do |v2, j|
					#puts "#{v1}:#{v2}"
					#puts hap_mat[i]
					#puts hap_mat[j]
					dist = DamerauLevenshtein.distance(hap_mat[i], hap_mat[j])
					#puts dist
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
end