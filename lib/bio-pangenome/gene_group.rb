module Bio
	module Pangenome
	class GeneGroup < Hash
		attr_accessor :gene 
		attr :sequences
		attr_accessor :expected_varieties



		def aligned_sequences
			return @aligned_sequences if @aligned_sequences
			return Hash.new if self.length == 0
			options = [ '--quiet', "--thread" , "8"]
			mafft = Bio::MAFFT.new( "mafft" , options)
			report = mafft.query_align(sequences)
			@aligned_sequences = report.alignment
			#@aligned_sequences.include AlignmentExtension
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



	end 

	class GeneGroupSet < Hash 
		attr_accessor :varieties

		def complete
			ret = GeneGroupSet.new { |h, k| h[k] = GeneGroup.new }
			self.each_pair do |k,v|
				ret[k] = v if v.keys.length == varieties.length
			end 
			return ret
		end
	end

	def self.load_sequences( varieties:[], genes:{}, prefix: "../flanking/filtered/",  suffix: ".cds.fa.gz", set_id: "cds" )
		ret = GeneGroupSet.new { |h, k| h[k] = GeneGroup.new }
		ret.varieties = varieties
		varieties.each do |variety|
			path = "#{prefix}/#{variety}#{suffix}"
			infile = open(path)
			io = Zlib::GzipReader.new(infile) 
			Bio::FlatFile.open(Bio::FastaFormat, io) do |fasta_file|
				fasta_file.each do |entry|
					definition, region = entry.definition.split("::")
					seq_name = parseSequenceName(region, definition) if region
					
					seq_name = GeneFlankingRegion.new(entry.definition,
						nil, "",
						"", entry.definition, set_id, nil, variety ) unless region
					arr = definition.split(".")
					next unless genes.include? seq_name.gene
					row = genes[seq_name.gene]
					row = genes[arr[0]] unless region
					seq_name.gene = row["gene"]
					
					seq = entry.seq
					seq.gsub!(/N*$/, '')
					seq.gsub!(/^N*/, '')
					seq_name.sequence = seq
					base_gene = row["gene"]
					ret[base_gene][variety] = seq_name unless ret[base_gene][variety]
					ret[base_gene].gene = base_gene
					ret[base_gene].expected_varieties = varieties
 				end
			end
			io.close
		end
		ret
	end

	end

end