module Bio::Pangenome
	class GeneGroup < Hash
		attr_accessor :gene 
		attr :sequences
		attr_accessor :expected_varieties

		def mafft 
			@@mafft ||= Bio::MAFFT.new( "mafft" , [ '--quiet', "--thread" , "4"])
		end

		def aligned_sequences
			return @aligned_sequences if @aligned_sequences
			return Hash.new if self.length == 0
			report = mafft.query_align(sequences)
			@aligned_sequences = report.alignment
			@aligned_sequences
    	end

    	def lengths
    		sequences.keys.map do |k|
    			sequences[k].length
    		end
    	end

    	def sequences
			return @sequences if @sequences 
			@sequences = Bio::Alignment::SequenceHash.new
		  	self.each_pair do |k,v|
		  		@sequences[k] = v.sequence	
		  	end
			@sequences
    	end

    	def mask 
    		@mask ||= Bio::Pangenome::HaplotypeMask.new(aligned_sequences)
    	end

    	def snps_hap
    		return @snps_hap if @snps_hap
    		@snps_hap = Hash.new 
    		expected_varieties.each do |variety|
    			@snps_hap[variety] = snps.map { |i| aligned_sequences[variety][i]}.join("")  
    		end
    		@snps_hap
    	end

    	def snps
    		mask.valid_snps
    	end

    	def hap_position(pos)
    		mask.valid_snps[pos]
    	end

    	def haplotype_matrix
    		matrix = []
    		snps.each do |s|
    			row  = [gene, s]
    			row.concat expected_varieties.map { |v| aligned_sequences[v][s] }
    			matrix << row
    		end
    		matrix
    	end


	end 

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
			hap.select{|h| not h.nil?  }.join(sep)
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