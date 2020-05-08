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

end